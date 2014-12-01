import os
import shutil
from jinja2 import Template

class Controller:
    def __init__(self, conf):
        self.conf = conf


    def validates(self):
        pass


    def _createJobDirs(self, jobdir, jobname, colvar, spec):
        inputdir = os.path.join(jobdir, 'input')
        outputdir = os.path.join(jobdir, 'output')
        colvarfile = os.path.join(inputdir, 'colvars.conf')

        if not os.path.exists(jobdir): os.makedirs(jobdir)
        if not os.path.exists(inputdir): os.makedirs(inputdir)
        if not os.path.exists(outputdir): os.makedirs(outputdir)

        colvar = self.colvars.append(colvar, jobname, **spec)
        colvar.write(colvarfile)

        shutil.copyfile(self.conf['psffile'], os.path.join(inputdir, self.conf['psffile']))
        shutil.copyfile(self.conf['pdbfile'], os.path.join(inputdir, self.conf['pdbfile']))

        for i in range(spec['nwin']):
            windir = os.path.join(outputdir, str(i))
            if not os.path.exists(windir): os.mkdir(windir)

        return jobdir, inputdir, outputdir

    def _createRMSD(self, jobname, jobdir, spec):
        jobdir, inputdir, outputdir = self._createJobDirs(jobdir, jobname, 'RMSD', spec)
        shutil.copyfile(spec['refpdb'], os.path.join(inputdir, spec['refpdb']))

    def _createAngle(self, jobname, jobdir, spec):
        jobdir, inputdir, outputdir = self._createJobDirs(jobdir, jobname, 'Angle', spec)

    def _createSepR(self, jobname, jobdir, spec):
        jobdir, inputdir, outputdir = self._createJobDirs(jobdir, jobname, 'Distance', spec)


    def create(self):
        for jobname in ['RMSDs/bb-receptor', 'RMSDs/bb-ligand', 'RMSDs/sc-receptor', 'RMSDs/sc-ligand',
                        'Angles/phi', 'Angles/theta', 'Angles/alpha', 'Angles/beta', 'Angles/gamma',
                        'Distance/sepR']:
            jobtype, jobname = jobname.split('/')
            basedir = os.path.join(self.conf['workdir'], jobtype)
            jobdir = os.path.join(basedir, jobname)
            spec = self.conf['simulations'][jobtype][jobname]

            if not spec.has_key('nwin'):
                spec['nwin'] = int((spec['max'] - spec['min']) / spec['bin']) + 1
            if not spec.has_key('numruns'):
                spec['numruns'] = self.conf['numruns']
            if not spec.has_key('temperature'):
                spec['temperature'] = self.conf['temperature']
            if not spec.has_key('jobname'):
                spec['jobname'] = jobname

            if jobtype == 'RMSDs':
                self._createRMSD(jobname, jobdir, spec)
            if jobtype == 'Angles':
                spec['centers'] = self.conf['simulations']['refvalues'][jobname]
                spec['refatoms'] = {}
                for ref in spec['angle']:
                    spec['refatoms'][ref] = self.conf['simulations']['refatoms'][ref]
                self._createAngle(jobname, jobdir, spec)
            if jobtype == 'Distance':
                spec['refatoms'] = {}
                for ref in spec['distance']:
                    spec['refatoms'][ref] = self.conf['simulations']['refatoms'][ref]
                self._createSepR(jobname, jobdir, spec)

            for filename in ('umbrella.namd', 'base.conf', 'sort.py', 'prepare.py'):
                realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'static', filename)
                shutil.copyfile(realname, os.path.join(jobdir, os.path.basename(realname)))

            for filename in ('remd.conf',):
                realname = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'templates', filename)
                filename = os.path.join(os.path.join(jobdir, os.path.basename(realname)))
                template = Template(open(realname).read())
                open(filename, 'w').write(template.render(spec))


    def validate(self):
        # validate selection
        pass