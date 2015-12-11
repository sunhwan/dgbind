import argparse
import yaml
import os, sys
from colvars import Colvars
import shutil
from jinja2 import Template
from MDAnalysis.coordinates.core import get_writer_for

def createJobDirs(jobdir, jobname, colvar, spec):
    inputdir = os.path.join(jobdir, 'input')
    outputdir = os.path.join(jobdir, 'output')
    colvarfile = os.path.join(inputdir, 'colvars.conf')

    if not os.path.exists(jobdir): os.makedirs(jobdir)
    if not os.path.exists(inputdir): os.makedirs(inputdir)
    if not os.path.exists(outputdir): os.makedirs(outputdir)

    spec['receptor'] = conf['receptor']['selection']
    spec['ligand'] = conf['ligand']['selection']
    colvar = colvars.append(colvar, jobname, **spec)
    colvar.write(colvarfile)

    for cv in colvars:
        if cv.spec.has_key('refpdb') and cv.spec.has_key('selected_atoms'):
            refpdb = os.path.join(inputdir, os.path.basename(cv.spec['refpdb']))
            writer = get_writer_for(format='PDB', multiframe=False)
            writer(refpdb).write(cv.spec['selected_atoms'])

    psffile = os.path.basename(spec['psffile'])
    pdbfile = os.path.basename(spec['pdbfile'])
    shutil.copyfile(spec['psffile'], os.path.join(inputdir, psffile))
    shutil.copyfile(spec['pdbfile'], os.path.join(inputdir, pdbfile))

    for i in range(spec['num_windows']):
        windir = os.path.join(outputdir, str(i))
        if not os.path.exists(windir): os.mkdir(windir)

    for filename in ('umbrella.namd', 'base.conf', 'sort.py', 'prepare.py'):
        realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'static', filename)
        shutil.copyfile(realname, os.path.join(jobdir, os.path.basename(realname)))

    spec['pbs'] = conf['pbs']
    for filename in ('remd.conf', 'namd.conf', 'run.pbs'):
        realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', filename)
        filename = os.path.join(os.path.join(jobdir, os.path.basename(realname)))
        template = Template(open(realname).read())
        open(filename, 'w').write(template.render(spec))

    return jobdir, inputdir, outputdir

def createRMSDs(jobname, jobdir, spec, colvars):
    if not spec.has_key('selection'):
        for k in ('bb-receptor', 'sc-receptor', 'bb-ligand', 'sc-ligand'):
            if jobname.startswith(k) and conf['refatoms'].has_key(k):
                spec['selection'] = conf['refatoms'][k]
    jobdir, inputdir, outputdir = createJobDirs(jobdir, jobname, 'RMSD', spec)
    if spec.has_key('refpdb'):
        refpdb = os.path.basename(spec['refpdb'])
        shutil.copyfile(spec['refpdb'], os.path.join(inputdir, refpdb))

def createAngles(jobname, jobdir, spec, colvars):
    spec['centers'] = spec['refvalue']
    spec['refatoms'] = {}
    for ref in spec['angle']:
        spec['refatoms'][ref] = conf['refatoms'][ref]
    jobdir, inputdir, outputdir = createJobDirs(jobdir, jobname, 'Angle', spec)

def createDistance(jobname, jobdir, spec, colvars):
    spec['refatoms'] = {}
    for ref in spec['distance']:
        spec['refatoms'][ref] = conf['refatoms'][ref]
    jobdir, inputdir, outputdir = createJobDirs(jobdir, jobname, 'Distance', spec)

def createSimulationInput(jobnames, colvars, psffile, pdbfile):
    for jobname in jobnames:
        jobtype, jobname = jobname.split('/')
        if not conf['simulations'].has_key(jobtype) or \
           not conf['simulations'][jobtype].has_key(jobname):
            continue
        basedir = os.path.join(conf['spec']['workdir'], jobtype)
        jobdir = os.path.join(basedir, jobname)
        spec = conf['simulations'][jobtype][jobname]

        for k in conf['spec'].keys():
            spec[k] = conf['spec'][k]
        for k in conf['simulations'][jobtype]['default'].keys():
            spec[k] = conf['simulations'][jobtype]['default'][k]
        spec['psffile'] = psffile
        spec['pdbfile'] = pdbfile

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

    default_jobs = ['RMSDs/bb-receptor-site',
                    'RMSDs/bb-ligand-site',
                    'RMSDs/sc-receptor-site',
                    'RMSDs/sc-ligand-site',
                    'Angles/phi',
                    'Angles/theta',
                    'Angles/alpha',
                    'Angles/beta',
                    'Angles/gamma',
                    'Distance/sepR']
    psffile = conf['complex']['psf']
    pdbfile = conf['complex']['pdb']
    colvars = Colvars(conf, psffile, pdbfile)
    colvars.append('Omega', 'Omega', selection='receptor and name C', receptor=conf.get('receptor').get('selection'))
    colvars.append('Pin', 'Pin', selection='receptor and name C', receptor=conf.get('receptor').get('selection'))
    createSimulationInput(default_jobs, colvars, psffile, pdbfile)

    psffile = conf['receptor']['psf']
    pdbfile = conf['receptor']['pdb']
    bulk_receptor_jobs = ['RMSDs/bb-receptor-bulk', 'RMSDs/sc-receptor-bulk']
    colvars = Colvars(conf, psffile, pdbfile)
    colvars.append('Omega', 'Omega', selection='receptor and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    colvars.append('Pin', 'Pin', selection='receptor and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    createSimulationInput(bulk_receptor_jobs, colvars, psffile, pdbfile)

    psffile = conf['ligand']['psf']
    pdbfile = conf['ligand']['pdb']
    bulk_ligand_jobs = ['RMSDs/bb-ligand-bulk', 'RMSDs/sc-ligand-bulk']
    colvars = Colvars(conf, psffile, pdbfile)
    colvars.append('Omega', 'Omega', selection='ligand and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    colvars.append('Pin', 'Pin', selection='ligand and name C', receptor=conf.get('receptor').get('selection'), ligand=conf.get('ligand').get('selection'))
    createSimulationInput(bulk_ligand_jobs, colvars, psffile, pdbfile)


