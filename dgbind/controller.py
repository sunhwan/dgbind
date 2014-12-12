from colvars import Colvars
import os

class Controller(object):
    def __init__(self, conf):
        self.conf = conf

    def _import_method(self, spec):
        method = self.conf['simulations']['method'] if not spec.has_key('method') else spec['method']
        name = "dgbind.methods.%s.controller" % method
        mod = __import__(name)
        components = name.split('.')
        for comp in components[1:]:
            mod = getattr(mod, comp)
        return mod

    def _createRMSDs(self, jobname, jobdir, spec, colvars):
        if not spec.has_key('selection'):
            for k in ('bb-receptor', 'sc-receptor', 'bb-ligand', 'sc-ligand'):
                if jobname.startswith(k) and self.conf['refatoms'].has_key(k):
                    spec['selection'] = self.conf['refatoms'][k]
        method = self._import_method(spec)
        controller = method.Controller(self.conf, colvars)
        controller.createRMSD(jobname, jobdir, spec)

    def _createAngles(self, jobname, jobdir, spec, colvars):
        method = self._import_method(spec)
        controller = method.Controller(self.conf, colvars)
        spec['centers'] = self.conf['refvalues'][jobname]
        spec['refatoms'] = {}
        for ref in spec['angle']:
            spec['refatoms'][ref] = self.conf['refatoms'][ref]
        controller.createAngle(jobname, jobdir, spec)

    def _createDistance(self, jobname, jobdir, spec, colvars):
        method = self._import_method(spec)
        controller = method.Controller(self.conf, colvars)
        spec['refatoms'] = {}
        for ref in spec['distance']:
            spec['refatoms'][ref] = self.conf['refatoms'][ref]
        controller.createDistance(jobname, jobdir, spec)

    def _createSimulationInput(self, jobnames, colvars, psffile, pdbfile):
        for jobname in jobnames:
            jobtype, jobname = jobname.split('/')
            if not self.conf['simulations'].has_key(jobtype) or \
               not self.conf['simulations'][jobtype].has_key(jobname):
                continue
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
            if self.conf.has_key('receptor'):
                spec['receptor'] = self.conf['receptor']
            if self.conf.has_key('ligand'):
                spec['ligand'] = self.conf['ligand']
            spec['psffile'] = psffile
            spec['pdbfile'] = pdbfile

            method = getattr(self, '_create%s' % jobtype)
            method(jobname, jobdir, spec, colvars)

    def create(self):
        if not self.conf['simulations'].has_key('method'):
            self.conf['simulations']['method'] = 'REMD'

        if self.conf.has_key('receptor'):
            for k,refatom in self.conf['refatoms'].items():
                self.conf['refatoms'][k] = self.conf['refatoms'][k].replace('receptor', self.conf['receptor'])
        if self.conf.has_key('ligand'):
            for k,refatom in self.conf['refatoms'].items():
                self.conf['refatoms'][k] = self.conf['refatoms'][k].replace('ligand', self.conf['ligand'])

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
        psffile = self.conf['psffile']['complex']
        pdbfile = self.conf['pdbfile']['complex']
        colvars = Colvars(self.conf, psffile, pdbfile)
        colvars.append('Omega', 'Omega', selection='receptor and name C', receptor=self.conf.get('receptor'), ligand=self.conf.get('ligand'))
        colvars.append('Pin', 'Pin', selection='receptor and name C', receptor=self.conf.get('receptor'), ligand=self.conf.get('ligand'))
        self._createSimulationInput(default_jobs, colvars, psffile, pdbfile)

        psffile = self.conf['psffile']['receptor']
        pdbfile = self.conf['pdbfile']['receptor']
        bulk_receptor_jobs = ['RMSDs/bb-receptor-bulk', 'RMSDs/sc-receptor-bulk']
        colvars = Colvars(self.conf, psffile, pdbfile)
        colvars.append('Omega', 'Omega', selection='receptor and name C', receptor=self.conf.get('receptor'), ligand=self.conf.get('ligand'))
        colvars.append('Pin', 'Pin', selection='receptor and name C', receptor=self.conf.get('receptor'), ligand=self.conf.get('ligand'))
        self._createSimulationInput(bulk_receptor_jobs, colvars, psffile, pdbfile)

        psffile = self.conf['psffile']['ligand']
        pdbfile = self.conf['pdbfile']['ligand']
        bulk_ligand_jobs = ['RMSDs/bb-ligand-bulk', 'RMSDs/sc-ligand-bulk']
        colvars = Colvars(self.conf, psffile, pdbfile)
        colvars.append('Omega', 'Omega', selection='ligand and name C', receptor=self.conf.get('receptor'), ligand=self.conf.get('ligand'))
        colvars.append('Pin', 'Pin', selection='ligand and name C', receptor=self.conf.get('receptor'), ligand=self.conf.get('ligand'))
        self._createSimulationInput(bulk_ligand_jobs, colvars, psffile, pdbfile)
