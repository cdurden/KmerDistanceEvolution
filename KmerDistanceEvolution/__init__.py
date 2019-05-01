# -*- coding: utf-8 -*-
"""

This module implements the XTree and the PhylogeneticTree classes used in the Coalescent,
and methods to test for isomorphisms between them XTree objects, to get the set of splits, and to test for refinement relations

Todo:
    * Finish documentation

.. _Google Python Style Guide:
   http://google.github.io/styleguide/pyguide.html

"""

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix
import numpy as np
import itertools

class XTree(Tree):
    """ An XTree implements X-tree, which is a Tree along with a labeling of its vertices.

    XTree is a subclass of Tree from Bio.Phylo.BaseTree. X-tree's are defined mathematically in Semple and Steel. Note: Tree from Bio.Phylo.BaseTre does not satisfy the requirement (in Semple and Steel) that vertices of degree <= 2 are labeled.

    The method used to test whether X-trees are isomorphic is to test for equality of their sets of splits. This method works even for a data structure with unlabeled vertices of degree <=2, because the unlabeled vertices do not affect splits.

    Properties created with the ``@property`` decorator should be documented
    in the property's getter method.

    Attributes:
        tree (:obj:`Tree`): A Tree.
        labeling (:obj:`dict`, optional): A dictionary representing the mapping between each Clade in tree and a label (str).
        X (:obj:`set`): The set of labels.
        total_order: For internal use only.

    """

    def __init__(self, tree, labeling=None):
        """ An XTree implements X-tree, which is a Tree along with a labeling of its vertices.
        """

        self.tree = tree
        if labeling is None:
            labeling = dict((clade,set([clade.name])) for clade in tree.get_terminals())
        assert(all([isinstance(labels,set) for labels in labeling.values()]))
        assert(all([labeled_clade in list(tree.find_clades()) for labeled_clade in labeling.keys()]))
        #assert(all([clade in labeling.keys() for clade in tree.find_clades() if degree(tree,clade)<=2]))
        self.labeling = labeling
        self.X = set(itertools.chain.from_iterable(labeling.values()))
        self.total_order = list()
        for clade in tree.find_clades():
            self.total_order += [clade]
    def get_path_index(self, fr, to):
        return([self.total_order.index(clade) for clade in fr.get_path(to)])
    def get_clade_labels(self, clade):
        labels = set()
        try:
            labels = labels.union(self.labeling[clade])
        except KeyError:
            pass
        for descendant in clade:
            labels = labels.union(self.get_clade_labels(descendant))
        return(labels)
    def get_splits(self, clade=None, include_trivial=False):
        splits = list()
        for clade in self.tree.find_clades():
            X0 = self.get_clade_labels(clade)
            if (len(X0) > 1 and len(self.X-X0) > 1) or include_trivial:
                if min(X0) < min(self.X-X0):
                    splits.append((tuple(sorted(X0)),tuple(sorted(self.X-X0))))
                else:
                    splits.append((tuple(sorted(self.X-X0)),tuple(sorted(X0))))
                    #splits.append((self.X-X0,X0))
#        if clade is None:
#            clade = self.tree.root
#        splits = list()
#        for descendant in clade:
#            X0 = self.get_clade_labels(descendant)
#            if len(X0) > 1 or include_trivial:
#                splits.append([X0,self.X-X0])
#            splits += self.get_splits(descendant)
        return(tuple(sorted(list(set(splits)))))
    def restricted_subtree(self, X):
        labeling = dict([(k,v.intersection(X)) for k,v in self.labeling.iteritems() if not v.isdisjoint(X)])
        restricted_subtree = XTree(self.tree,labeling)
        return(restricted_subtree)
    def refines(self, xtree):
        splits_of_self = self.get_splits()
        splits_of_xtree = xtree.get_splits()
        for split_of_self in splits_of_self:
            # Splits have two blocks, so we test if one of the blocks (the first one) of split_of_self matches either block of any split_of_xtree
            if not any([split_of_self[0] in split_of_xtree for split_of_xtree in splits_of_xtree]):
                return(False)
        return(True)
    def displays(self, xtree):
        # Check whether self displays xtree
        # That is, check whether xtree refines the restriction of self to the label set of xtree (xtree.X)
        return(xtree.refines(self.restricted_subtree(xtree.X)))
    def is_isomorphic(self, xtree):
        return(self.refines(xtree) and xtree.refines(self))
    

class PhylogeneticTree(XTree):
    def __init__(self, tree, labeling):
        self.tree = tree
        assert(all([isinstance(labels,set) for labels in labeling.values()]))
        assert(all([labeled_clade in list(tree.find_clades()) for labeled_clade in labeling.keys()]))
        for clade in tree.find_clades():
            nr_children = len(clade.clades)
            assert(nr_children==0 or nr_children==2)
        self.labeling = labeling
        self.X = set(itertools.chain.from_iterable(labeling.values()))
        self.total_order = list()
        for clade in tree.find_clades():
            self.total_order += [clade]

