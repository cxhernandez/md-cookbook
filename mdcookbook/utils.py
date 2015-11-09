from mdcookbook import __author__

from simtk.openmm import XmlSerializer
import argparse
import time
import os


class timing(object):
    "Context manager for printing performance"

    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = time.time()

    def __exit__(self, ty, val, tb):
        end = time.time()
        print("PERFORMANCE [%s] : %0.3f seconds" % (self.name,
                                                    end - self.start))
        return False


def serialize(obj, dirname, objname):
    filename = './%s/%s' % (dirname, objname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)
    with open(filename, 'wb') as objfile:
        objfile.write(XmlSerializer.serialize(obj))


def count(obj):
    for i, _ in enumerate(obj):
        pass
    return i


def get_args():
    parser = argparse.ArgumentParser(
        epilog="Written by %s" % __author__,
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '-pt',
        '--platform',
        dest='platform',
        help='Platform type.',
        choices=[
            'CPU',
            'CUDA',
            'OpenCL'],
        default="CPU")
    parser.add_argument('-d', '--device-index', dest='device',
                        help='GPU device index.', default=0, type=int)
    parser.add_argument('-nc', '--n-clones', dest='n_clones',
                        help='Number of clones to create.', default=5,
                        type=int)
    parser.add_argument('-n', '--n-steps', dest='n_steps',
                        help='Number of equilibration steps.',
                        default=50000, type=int)
    parser.add_argument('-i', '--max-iterations', dest='max_iter',
                        help='Max iterations for minimization.',
                        default=1000, type=int)
    parser.add_argument('-c', '--ion-concentation', dest='ion_content',
                        help='Ion concentration (molar).',
                        default=0.1, type=float)
    parser.add_argument('-t', '--temperature', dest='temp',
                        help='Simulation temperature (kelvin).',
                        default=300, type=int)
    parser.add_argument('-b', '--box-size', dest='box_size',
                        help='Boxsize expressed as: Vec3(x, y, z).',
                        default='Vec3(5, 5, 5)')
    parser.add_argument('-w', '--n-waters', dest='n_wat',
                        help='Number of waters in solvent '
                        '(Overrides BOX_SIZE).',
                        default=None, type=int)
    return parser
