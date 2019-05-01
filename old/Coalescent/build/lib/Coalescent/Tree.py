from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix
import numpy as np
import itertools

class XTree(Tree):
    # Note: An XTree is represented by a tree data structure along with a labeling of its vertices.
    # The data structure does not satisfy the requirement that vertices of degree <= 2 are labeled.
    # We get a real X-tree by suppressing unlabeled vertices of degree <= 2.
    # You might say that this implementation is based on the definition of an X-tree as an equivalence class of trees where trees are equivalent if the trees obtained by suppressing all unlabeled vertices of degree <= 2 are isomorphic.
    # This representation is okay for the purposes of testing whether X-trees are isomorphic because the unlabeled vertices do not affect splits.
    def __init__(self, tree, labeling=None):
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

def huelsenbeck_tree(t_a,t_b,n):
    import numpy as np
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

