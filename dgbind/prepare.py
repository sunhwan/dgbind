import argparse
import yaml
import os, sys

def validate_configuration(conf):
    keys = ['receptor', 'ligand', 'psffile', 'pdbfile', 'workdir']
    msg = "%s section is missing in the configuration file."
    for k in keys:
        if not conf.has_key(k): return False, msg % k
    return True, ''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('conf', metavar='conf', help='Configuration file')
    #parser.add_argument('psf', metavar='psffile', help='PSF file')
    #parser.add_argument('pdb', metavar='pdbfile', help='PDB file')
    #parser.add_argument('-o', dest='output', help='output CRD file', default=None)
    args = parser.parse_args()

    if not os.path.exists(args.conf):
        print "Configuration file not found"
        sys.exit()

    conf = yaml.load(open(args.conf).read())
    is_valid, msg = validate_configuration(conf)
    if not is_valid:
        print "Configuration file has an error:"
        print msg
        sys.exit()

    from methods.controller import Controller
    controller = Controller(conf)
    controller.create()
