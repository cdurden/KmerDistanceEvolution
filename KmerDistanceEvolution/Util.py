# -*- coding: utf-8 -*-
"""Simulation of coalescence and mutations

.. moduleauthor:: Chris Durden <cdurden@gmail.com>

"""
import shutil
import numpy as np
import tempfile
from . import Tree, Clade

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
    """ Assigns (arbitrary) names to clades of tree """
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

def huelsenbeck_tree(t_a,t_b,n):
    import numpy as np
    from Bio import Phylo
    from cStringIO import StringIO
    if n==4:
        labels = ['A','B','C','D']
        np.random.shuffle(labels)
        #treedata = "((A:{0:f}, B:{1:f}):{0:f}, (C:{0:f},D:{1:f}):0):0".format(t_a,t_b)
        treedata = "(({2:s}:{0:f}, {3:s}:{1:f}):{0:f}, ({4:s}:{0:f},{5:s}:{1:f}):0):0".format(t_a,t_b,*labels)
        handle = StringIO(treedata)
        tree = Phylo.read(handle, "newick")
        return(tree)
    elif n==5:
        labels = ['A','B','C','D','E']
        np.random.shuffle(labels)
        #treedata = "(((A:{0:f}, B:{1:f}):{0:f},C:{1:f}):0, (D:{0:f},E:{1:f}):{0:f}):0".format(t_a,t_b)
        treedata = "((({2:s}:{0:f}, {3:s}:{1:f}):{0:f},{4:s}:{1:f}):0, ({5:s}:{0:f},{6:s}:{1:f}):{0:f}):0".format(t_a,t_b,*labels)
        handle = StringIO(treedata)
        tree = Phylo.read(handle, "newick")
        return(tree)

def enumerate_5_taxon_tree_labelings(labels=None):
    if labels is None:
        labels = set([0,1,2,3,4])
    else:
        labels = set(tuple(labels))
        assert(len(labels)==5)
    labelings = list()
    for i in labels:
        for j in labels-set([i,min(labels-{i})]):
            #labelings.append([i,min(labels-{i}),j])
            labelings.append([min(labels-{i}),j,i])
            labelings[-1] += sorted(labels-set(labelings[-1]))
    labelings = [tuple(labeling) for labeling in labelings]
    return(labelings)

def enumerate_5_taxon_trees(s_terminals):
    assert(len(s_terminals)==5)
    labelings = enumerate_5_taxon_tree_labelings()
    trees = list()
    for labeling in labelings:
         yield make_tree([s_terminals[i] for i in labeling])

def make_tree(labels):
    if len(labels)==1:
        return(Tree.from_clade(Clade(name=labels[0])))
    else:
        return(Tree.from_clade(Clade(clades=[make_tree(labels[:-1]).root,Clade(name=labels[-1])])))

