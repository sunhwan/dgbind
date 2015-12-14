from flask import Flask, request, render_template, redirect, url_for, send_file
from werkzeug import secure_filename
import os
import yaml

app = Flask(__name__)
app.config['upload_folder'] = '/tmp'
app.config['job_folder'] = '/tmp'
app.debug = True

def get_job_folder(uuid):
    return os.path.join(app.config['job_folder'], uuid)

def append_ref_atoms(psffile, pdbfile, reffile, receptor_segid, ligand_segid):
    import MDAnalysis as mda
    import StringIO
    u = mda.Universe(psffile, pdbfile)
    u.selectAtoms('segid %s or segid %s' % (receptor_segid, ligand_segid)).write(reffile, bonds=None)

    ref = StringIO.StringIO()
    ref.writelines(open(reffile).readlines()[:-1])
    ref.write("""ATOM      1  DUM P1  Z   1       0.000   0.000   0.000  1.00  0.00      DUM  C 0
ATOM      2  DUM P2  Z   2       0.000   0.000   0.000  1.00  0.00      DUM  C 0
ATOM      3  DUM P3  Z   3       0.000   0.000   0.000  1.00  0.00      DUM  C 0
ATOM      4  DUM L1  Z   4       0.000   0.000   0.000  1.00  0.00      DUM  C 0
ATOM      5  DUM L2  Z   5       0.000   0.000   0.000  1.00  0.00      DUM  C 0
ATOM      6  DUM L3  Z   6       0.000   0.000   0.000  1.00  0.00      DUM  C 0""")
    ref.seek(0)
    open(reffile, 'w').write(ref.read())

def pick_reference_atoms(psffile, pdbfile, receptor_segid, ligand_segid, job_basedir, selections=[]):
    reffile = os.path.join(job_basedir, 'reference.pdb')
    if not os.path.exists(reffile) or 1:
        # append reference atoms to the target pdbfile
        append_ref_atoms(psffile, pdbfile, reffile, receptor_segid, ligand_segid)

    import MDAnalysis as mda
    u = mda.Universe(reffile)
    if not selections:
        import random
        # P1
        selections.append('segid %s and name CA' % receptor_segid)

        # P2
        first_rid = u.atoms.selectAtoms('segid %s' % receptor_segid)[0].resid
        last_rid = u.atoms.selectAtoms('segid %s' % receptor_segid)[-1].resid
        rid = random.randrange(first_rid, last_rid)
        selections.append('segid %s and resid %d and name CA' % (receptor_segid, rid))

        # P3
        rid = random.randrange(first_rid, last_rid)
        selections.append('segid %s and resid %d and name CA' % (receptor_segid, rid))

        # L1
        selections.append('segid %s and name CA' % ligand_segid)

        # L2
        first_rid = u.atoms.selectAtoms('segid %s' % ligand_segid)[0].resid
        last_rid = u.atoms.selectAtoms('segid %s' % ligand_segid)[-1].resid
        rid = random.randrange(first_rid, last_rid)
        selections.append('segid %s and resid %d and name CA' % (ligand_segid, rid))

        # L3
        rid = random.randrange(first_rid, last_rid)
        selections.append('segid %s and resid %d and name CA' % (ligand_segid, rid))

    ref_atoms = []
    for i in range(6):
        ref_atoms.append(u.atoms.selectAtoms('segid DUM and resid %d' % (i+1)))
        ref_atoms[i].set_positions(u.atoms.selectAtoms(selections[i]).centerOfMass())
    u.atoms.write(reffile, bonds=None)

    alpha = (ref_atoms[1] + ref_atoms[0] + ref_atoms[3]).angle()
    beta  = (ref_atoms[2] + ref_atoms[1] + ref_atoms[0] + ref_atoms[3]).dihedral()
    gamma = (ref_atoms[1] + ref_atoms[0] + ref_atoms[3] + ref_atoms[4]).dihedral()
    r     = (ref_atoms[0] + ref_atoms[3]).bond()
    theta = (ref_atoms[0] + ref_atoms[3] + ref_atoms[4]).angle()
    phi   = (ref_atoms[0] + ref_atoms[3] + ref_atoms[4] + ref_atoms[5]).dihedral()

    return selections, {'alpha': alpha, 'beta': beta, 'gamma': gamma, 'r': r, 'theta': theta, 'phi': phi}

def read_cached_conf(jobid, keys=[]):
    job_basedir = get_job_folder(jobid)
    conf_filename = os.path.join(job_basedir, 'conf.yaml')
    if os.path.exists(conf_filename):
        conf = yaml.load(open(os.path.join(job_basedir, 'conf.yaml')).read())
    else:
        conf = {}

    for k in keys:
        conf[k] = request.form[k]
    open(conf_filename, 'w').write(yaml.dump(conf))
    return conf

def prepare_input_files(jobid, conf):
    job_basedir = get_job_folder(jobid)


@app.route('/download/<jobid>/<filename>')
@app.route('/download/<jobid>', defaults={'filename': None})
def download(jobid, filename):
    job_basedir = get_job_folder(jobid)
    if not os.path.exists(job_basedir): abort(404)

    if filename:
        filename = os.path.join(job_basedir, os.path.basename(filename))
        if not os.path.exists(filename): abort(404)
        fp = open(filename)
        if filename.endswith('gif'):
            return send_file(fp, as_attachment=False, attachment_filename=os.path.basename(filename))
        return send_file(fp, as_attachment=True, attachment_filename=os.path.basename(filename))

    else:
        import tarfile
        import StringIO, glob
        fp = StringIO.StringIO()
        tar = tarfile.open(fileobj=fp, mode='w:gz')
        excludes = ['taskmanager.py', 'err', 'out', 'run.pbs']
        for f in glob.glob('%s/*' % job_basedir):
            if os.path.basename(f) in excludes: continue
            tar.add(f, arcname='ppiserver/%s' % os.path.basename(f))
        tar.close()
        fp.seek(0)
        return send_file(fp, mimetype='application/x-gzip', as_attachment=True, attachment_filename='ppiserver.tar.gz')


@app.route("/input-files/<jobid>", methods=["POST"])
def generate(jobid):
    # user defined simulation setting and cache settings
    conf = read_cached_conf(jobid, ('rmsd[receptor][selection]', 'rmsd[ligand][selection]', 'rmsd[method]', 'rmsd[nwin]', 'rmsd[bin]',
                                    'sc[receptor][selection]', 'sc[ligand][selection]', 'sc[method]', 'sc[nwin]', 'sc[bin]',
                                    'angles[method]', 'angles[nwin]', 'angles[bin]',
                                    'sepr[method]', 'sepr[nwin]', 'sepr[bin]', 'sepr[rmin]'))
    prepare_input_files(jobid, conf)

    return render_template('inputs.html', jobid=jobid)


@app.route("/define-method/<jobid>", methods=["POST"])
def method(jobid):
    # user defined simulation setting and cache settings
    conf = read_cached_conf(jobid, ('selections[0]', 'selections[1]', 'selections[2]', 'selections[3]', 'selections[4]', 'selections[5]'))

    methods_available = (('umbrella', 'Umbrella Sampling'),
                         ('umbrella-remd', 'Umbrella Sampling with Replica Exchange'),
                         ('abf', 'Adaptive Biasing Force'), )
    default_method = 'umbrella-remd'

    return render_template('method.html', jobid=jobid, methods_available=methods_available, default_method=default_method)
    

@app.route("/define-reference-atoms/<jobid>", methods=["POST"])
def reference(jobid):
    # user defined simulation setting and cache settings
    conf = read_cached_conf(jobid, ('psffile', 'pdbfile', 'jobid', 'temperature', 'pressure', 'fftype', 'receptor-segid', 'ligand-segid'))

    # suggest initial reference atoms
    job_basedir = get_job_folder(jobid)
    psffile = os.path.join(job_basedir, conf['psffile'])
    pdbfile = os.path.join(job_basedir, conf['pdbfile'])
    selections, refvalues = pick_reference_atoms(psffile, pdbfile, conf['receptor-segid'], conf['ligand-segid'], job_basedir)

    return render_template('reference.html', selections=selections, refvalues=refvalues, jobid=jobid)


@app.route("/upload", methods=['GET', "POST"])
def upload():
    if request.method == 'GET': return redirect(url_for('index'))

    psf = request.files['psffile']
    pdb = request.files['pdbfile']

    if not psf and not pdb:
        error={'code': 500, 'message': 'There was problem uploading PDB files'}
        return render_template('index.html', error=error)

    psffile = secure_filename(psf.filename)
    pdbfile = secure_filename(pdb.filename)
    psf.save(os.path.join(app.config['upload_folder'], psffile))
    pdb.save(os.path.join(app.config['upload_folder'], pdbfile))

    # TODO:
    # 1. test number of atom mismatch
    # 2. test if more than two segments exist
    segments = []
    def test_psf_pdb(psffile, pdbfile):
        import MDAnalysis
        u = MDAnalysis.Universe(psffile, pdbfile)
        for segment in u.segments:
            segments.append(segment.id)

    test_psf_pdb(os.path.join(app.config['upload_folder'], psffile), 
                 os.path.join(app.config['upload_folder'], pdbfile))

    # generate temporary disk space and relocate psf/pdb files to the directory
    import uuid
    jobid = str(uuid.uuid4())
    job_basedir = get_job_folder(jobid)
    os.mkdir(job_basedir)

    from shutil import copy
    copy(os.path.join(app.config['upload_folder'], psffile), os.path.join(job_basedir, psffile))
    copy(os.path.join(app.config['upload_folder'], pdbfile), os.path.join(job_basedir, pdbfile))

    return render_template('upload.html', pdbfile=pdbfile, psffile=psffile, segments=segments, jobid=jobid)


@app.route("/")
def index():
    return render_template('index.html')

if __name__ == "__main__":
    app.run()
