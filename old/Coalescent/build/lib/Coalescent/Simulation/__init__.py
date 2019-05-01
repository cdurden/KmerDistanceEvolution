import shutil
import numpy as np
import tempfile

class TemporaryDirectory(object):
    """Context manager for tempfile.mkdtemp() so it's usable with "with" statement."""
    def __init__(self, persist=False):
        self.persist = persist
    def __enter__(self):
        self.name = tempfile.mkdtemp()
        return self.name

    def __exit__(self, exc_type, exc_value, traceback):
        if not self.persist:
            shutil.rmtree(self.name)

def name_clades(tree):
    """ Assigns names to clades of tree """
    names = {}
    existing_names = [clade.name for clade in tree.find_clades()]
    numbers = [str(i) for i in range(len(list(tree.find_clades())))]
    unused_numbers = list(set(numbers)-set(existing_names))
    idx = 0
    for clade in tree.find_clades():
        if not clade.name:
            clade.name = "v{:s}".format(unused_numbers[idx])
            idx += 1
        names[clade.name] = clade
    return names
