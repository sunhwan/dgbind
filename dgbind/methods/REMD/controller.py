import os
import shutil
from jinja2 import Template
from MDAnalysis.coordinates.core import get_writer_for

class Controller(object):
    def __init__(self, conf, colvars):
        self.conf = conf
        self.colvars = colvars

    def validates(self):
        """Validates selection and input psf/pdb files"""
        pass

    def createJobDirs(self, jobdir, jobname, colvar, spec):
        inputdir = os.path.join(jobdir, 'input')
        outputdir = os.path.join(jobdir, 'output')
        colvarfile = os.path.join(inputdir, 'colvars.conf')

        if not os.path.exists(jobdir): os.makedirs(jobdir)
        if not os.path.exists(inputdir): os.makedirs(inputdir)
        if not os.path.exists(outputdir): os.makedirs(outputdir)

        colvar = self.colvars.append(colvar, jobname, **spec)
        colvar.write(colvarfile)

        for cv in self.colvars:
            if cv.spec.has_key('refpdb') and cv.spec.has_key('selected_atoms'):
                refpdb = os.path.join(inputdir, os.path.basename(cv.spec['refpdb']))
                writer = get_writer_for(format='PDB', multiframe=False)
                writer(refpdb).write(cv.spec['selected_atoms'])

        psffile = os.path.basename(spec['psffile'])
        pdbfile = os.path.basename(spec['pdbfile'])
        shutil.copyfile(spec['psffile'], os.path.join(inputdir, psffile))
        shutil.copyfile(spec['pdbfile'], os.path.join(inputdir, pdbfile))

        for i in range(spec['nwin']):
            windir = os.path.join(outputdir, str(i))
            if not os.path.exists(windir): os.mkdir(windir)

        for filename in ('umbrella.namd', 'base.conf', 'sort.py', 'prepare.py'):
            realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'static', filename)
            shutil.copyfile(realname, os.path.join(jobdir, os.path.basename(realname)))

        for filename in ('remd.conf',):
            realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', filename)
            filename = os.path.join(os.path.join(jobdir, os.path.basename(realname)))
            template = Template(open(realname).read())
            open(filename, 'w').write(template.render(spec))

        return jobdir, inputdir, outputdir

    def createRMSD(self, jobname, jobdir, spec):
        jobdir, inputdir, outputdir = self.createJobDirs(jobdir, jobname, 'RMSD', spec)
        refpdb = os.path.basename(spec['refpdb'])
        shutil.copyfile(spec['refpdb'], os.path.join(inputdir, refpdb))

    def createAngle(self, jobname, jobdir, spec):
        jobdir, inputdir, outputdir = self.createJobDirs(jobdir, jobname, 'Angle', spec)

    def createDistance(self, jobname, jobdir, spec):
        jobdir, inputdir, outputdir = self.createJobDirs(jobdir, jobname, 'Distance', spec)

