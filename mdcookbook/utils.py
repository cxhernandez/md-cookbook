from simtk.openmm import XmlSerializer
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
