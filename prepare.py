import argparse
import yaml
import os, sys
from dgbind.controller import Controller

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('conf', metavar='conf', help='Configuration file')
    args = parser.parse_args()

    if not os.path.exists(args.conf):
        print "Configuration file not found"
        sys.exit()

    conf = yaml.load(open(args.conf).read())
    #is_valid, msg = validate_configuration(conf)
    #if not is_valid:
    #    print "Configuration file has an error:"
    #    print msg
    #    sys.exit()

    controller = Controller(conf)
    controller.create()
