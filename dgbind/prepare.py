import argparse
import yaml
import os, sys
from colvars import Colvars
import shutil
from jinja2 import Template
import MDAnalysis as mda
import numpy as np
from MDAnalysis.coordinates.core import get_writer_for
import copy

def createJobDirs(jobdir, jobname, colvars, colvar_type, spec):
    inputdir = spec['inputdir']
    outputdir = spec['outputdir']
    colvarfile = spec['colvarfile']

    spec['receptor'] = conf['receptor']['selection']
    spec['ligand'] = conf['ligand']['selection']
    spec['jobname'] = jobname

    colvar = colvars.append(colvar_type, jobname, **spec)
    psffile = os.path.basename(spec['psffile'])
    pdbfile = os.path.basename(spec['pdbfile'])
    shutil.copyfile(spec['psffile'], os.path.join(inputdir, psffile))
    shutil.copyfile(spec['pdbfile'], os.path.join(inputdir, pdbfile))

    for prm in spec['parameters']:
        prmfile = os.path.basename(prm)
        shutil.copyfile(prm, os.path.join(inputdir, prmfile))

    spec['pbs'] = conf['pbs']

    # PME size
    universe = mda.Universe(spec['psffile'], spec['pdbfile'])
    bbox = universe.atoms.bbox()
    spec['A'], spec['B'], spec['C'] = bbox[1,:] - bbox[0,:]
    spec['xcen'], spec['ycen'], spec['zcen'] = np.average(bbox, axis=0)

    if spec.has_key('rest_selection'):
        sptpdb = os.path.join(inputdir, os.path.basename('spt.pdb'))
        writer = get_writer_for(format='PDB', multiframe=False)
        selected_atoms = universe.select_atoms(spec['rest_selection'])
        selected_atoms.set_bfactors(1)
        writer(sptpdb, bonds=False).write(universe.atoms)

    # initial structure preparation
    spec['use_rest'] = True
    for fname in ('input/namd.conf', 'input/run.pbs', 'input/prepare.py'):
        realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', fname)
        filename = os.path.join(os.path.join(inputdir, os.path.basename(realname)))
        if fname == 'input/run.pbs' and spec['jobtype'] == 'Angles':
            # angular restraints has minimum at the middle of the arranged windows
            realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', 'input/run-Angles.pbs')
        template = Template(open(realname).read())
        open(filename, 'w').write(template.render(spec, psffile=psffile, pdbfile=pdbfile))

    # input files preparation
    for i in range(spec['num_windows']):
        windir = os.path.join(outputdir, str(i))
        if not os.path.exists(windir): os.mkdir(windir)

    for filename in ('umbrella.namd', 'base.conf', 'sort.py', 'prepare.py'):
        realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'static', filename)
        shutil.copyfile(realname, os.path.join(jobdir, os.path.basename(realname)))

    spec['use_rest'] = False
    for filename in ('remd.conf', 'namd.conf', 'run.pbs'):
        realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', filename)
        filename = os.path.join(os.path.join(jobdir, os.path.basename(realname)))
        template = Template(open(realname).read())
        open(filename, 'w').write(template.render(spec))

    return jobdir, inputdir, outputdir, colvars

def createRMSDs(jobname, jobdir, spec, colvars):
    jobdir, inputdir, outputdir, colvars = createJobDirs(jobdir, jobname, colvars, 'RMSD', spec)
    if spec.has_key('refpdb'):
        refpdb = os.path.basename(spec['refpdb'])
        shutil.copyfile(spec['refpdb'], os.path.join(inputdir, refpdb))
    return colvars

def createAngles(jobname, jobdir, spec, colvars):
    spec['centers'] = spec['refvalue']
    spec['refatoms'] = {}
    for ref in spec['angle']:
        spec['refatoms'][ref] = conf['refatoms'][ref]
    jobdir, inputdir, outputdir, colvars = createJobDirs(jobdir, jobname, colvars, 'Angle', spec)
    return colvars

def createDistance(jobname, jobdir, spec, colvars):
    spec['centers'] = spec['refvalue']
    spec['refatoms'] = {}
    for ref in spec['distance']:
        spec['refatoms'][ref] = conf['refatoms'][ref]
    jobdir, inputdir, outputdir, colvars = createJobDirs(jobdir, jobname, colvars, 'Distance', spec)
    return colvars

def createSimulationInput(jobnames, colvars, psffile, pdbfile):
    for jobname in jobnames:
        print jobname
        jobtype, jobname = jobname.split('/')
        if not conf['simulations'].has_key(jobtype) or \
           not conf['simulations'][jobtype].has_key(jobname):
            continue
        basedir = os.path.join(conf['spec']['workdir'], jobtype)
        jobdir = os.path.join(basedir, jobname)
        inputdir = os.path.join(jobdir, 'input')
        outputdir = os.path.join(jobdir, 'output')
        colvarfile = os.path.join(inputdir, 'colvars.conf')

        if not os.path.exists(jobdir): os.makedirs(jobdir)
        if not os.path.exists(inputdir): os.makedirs(inputdir)
        if not os.path.exists(outputdir): os.makedirs(outputdir)

        spec = conf['simulations'][jobtype][jobname]
        spec['inputdir'] = inputdir
        spec['outputdir'] = outputdir
        spec['colvarfile'] = colvarfile
        spec['jobtype'] = jobtype
        spec['jobname'] = jobname
        spec['psffile'] = psffile
        spec['pdbfile'] = pdbfile

        for k in conf['spec'].keys():
            spec[k] = conf['spec'][k]
        for k in conf['simulations'][jobtype]['default'].keys():
            if not spec.has_key(k):
                spec[k] = conf['simulations'][jobtype]['default'][k]
            if spec.has_key('windows'):
                for rid in spec['windows'].keys():
                    if spec['windows'][rid].has_key(k): continue
                    spec['windows'][rid][k] = conf['simulations'][jobtype]['default'][k]
        for k in conf['namd'].keys():
            spec[k] = conf['namd'][k]
        if conf['receptor'].has_key('selection'):
            for k in ('selection', 'rest_selection'):
                if not spec.has_key(k): continue
                spec[k] = spec[k].replace('receptor', conf['receptor']['selection'])
        if conf['ligand'].has_key('selection'):
            for k in ('selection', 'rest_selection'):
                if not spec.has_key(k): continue
                spec[k] = spec[k].replace('ligand', conf['ligand']['selection'])

        if jobtype == 'Angles':
            angles = {'alpha': ['p1', 'l1', 'l2'],
                      'beta':  ['p1', 'l1', 'l2', 'l3'],
                      'gamma': ['p2', 'p1', 'l1', 'l2'],
                      'theta': ['p2', 'p1', 'l1'],
                      'phi': ['p3', 'p2', 'p1', 'l1']}
            spec['angle'] = angles[jobname]

        if jobtype == 'Distance':
            spec['distance'] = ['p1', 'l1']

        #method = getattr('create%s' % jobtype)
        method = globals()['create%s' % jobtype]
        method(jobname, jobdir, spec, colvars)

        colvars.write(colvarfile, inputdir=inputdir)
        for cv in colvars:
            if cv.spec.has_key('refpdb') and cv.spec.has_key('selected_atoms'):
                refpdb = os.path.join(inputdir, os.path.basename(cv.spec['refpdb']))
                writer = get_writer_for(format='PDB', multiframe=False)
                writer(refpdb).write(cv.spec['selected_atoms'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('conf', metavar='conf', help='Configuration file')
    args = parser.parse_args()

    if not os.path.exists(args.conf):
        print "Configuration file not found"
        sys.exit()

    conf = yaml.load(open(args.conf).read())

    if conf.has_key('receptor'):
        for k,refatom in conf['refatoms'].items():
            conf['refatoms'][k] = conf['refatoms'][k].replace('receptor', conf['receptor']['selection'])
    if conf.has_key('ligand'):
        for k,refatom in conf['refatoms'].items():
            conf['refatoms'][k] = conf['refatoms'][k].replace('ligand', conf['ligand']['selection'])

    default_jobs = ['RMSDs/bb_receptor_site',
                    'RMSDs/bb_ligand_site',
                    'RMSDs/sc_receptor_site',
                    'RMSDs/sc_ligand_site',
                    'Angles/alpha',
                    'Angles/beta',
                    'Angles/gamma',
                    'Angles/theta',
                    'Angles/phi',
                    'Distance/sepR']
    psffile = conf['complex']['psf']
    pdbfile = conf['complex']['pdb']
    colvars = Colvars(conf, psffile, pdbfile)
    colvars.append('Omega', 'Omega', selection='receptor and name C', receptor=conf.get('receptor').get('selection'))
    colvars.append('Pin', 'Pin', selection='receptor and name C', receptor=conf.get('receptor').get('selection'))
    createSimulationInput(default_jobs, colvars, psffile, pdbfile)

    psffile = conf['receptor']['psf']
    pdbfile = conf['receptor']['pdb']
    bulk_receptor_jobs = ['RMSDs/bb_receptor_bulk', 'RMSDs/sc_receptor_bulk']
    colvars = Colvars(conf, psffile, pdbfile)
    colvars.append('Omega', 'Omega', selection='receptor and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    colvars.append('Pin', 'Pin', selection='receptor and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    createSimulationInput(bulk_receptor_jobs, colvars, psffile, pdbfile)

    psffile = conf['ligand']['psf']
    pdbfile = conf['ligand']['pdb']
    bulk_ligand_jobs = ['RMSDs/bb_ligand_bulk', 'RMSDs/sc_ligand_bulk']
    colvars = Colvars(conf, psffile, pdbfile)
    colvars.append('Omega', 'Omega', selection='ligand and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    colvars.append('Pin', 'Pin', selection='ligand and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    createSimulationInput(bulk_ligand_jobs, colvars, psffile, pdbfile)


