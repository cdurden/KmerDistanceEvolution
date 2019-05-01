from .. import name_clades
from Bio.Phylo.BaseTree import Clade, Tree
from copy import copy
from itertools import combinations
import numpy as np

class ExponentialDistribution(object):
    def __init__(self, rate):
        self.rate = rate 
    def sample(self):
        dt = np.random.exponential(1.0/self.rate)
        return(dt)

class GeneLineage(Clade):
    pass

class EmbeddedGeneForest(object):
    # This is the object that coalesces to generate a gene tree (embedded within the species tree)
    def __init__(self, tree, embedding):
        self.tree = tree
        self.s_terminals = tree.get_terminals()
        assert(all([k in list(tree.find_clades()) for k in embedding.keys()]))
        self.embedding = dict((k,copy(v)) for k,v in embedding.items())
        for s_clade in tree.find_clades(order='postorder'):
            if not self.embedding.has_key(s_clade):
                self.embedding[s_clade] = []
        assert(self.is_sorted())
    def is_sorted(self):
        for s_clade in self.tree.find_clades():
            for i,j in combinations(range(len(self.embedding[s_clade])),2):
                if self.embedding[s_clade][i].is_parent_of(self.embedding[s_clade][j]) and i > j:
                    return False
        return True
    def maximal_gene_lineages(self, s_clade=None):
        # A gene lineage is a node in the gene forest
        # Given an s_clade (default: root of the tree), return the maximal gene_lineages which are embedded at or below (i.e. relative to the inverse partial order on the tree, where the root is the maximal element) the given s_clade
        if s_clade is None:
            s_clade = self.tree.root
        maximal_gene_lineages = set()
        for s_subclade in s_clade.find_clades(order='preorder'):
            for gene_lineage in self.embedding[s_subclade]:
                # We need the lists of clades that comprise the embedding to be sorted so that gene_lineages which are maximal, among those embedded within a given s_clade, are reached first. This should be the case because in the join method, we add the clade representing the coalescence to the embedding at index 0, before each of the coalescing clades.
                if not any([bool(maximal_gene_lineage.is_parent_of(gene_lineage)) for maximal_gene_lineage in maximal_gene_lineages]):
                    maximal_gene_lineages.add(gene_lineage)
        return(maximal_gene_lineages)
    def where(self, gene_lineage):
        for i,gene_lineages in enumerate(self.embedding.values()):
            if gene_lineage in gene_lineages:
                return(self.embedding.keys()[i])
    def join(self, clades, where, branch_length=None):
        if branch_length is not None:
            coalescence = Clade(clades=clades, branch_length=branch_length)
        else:
            coalescence = Clade(clades=clades)
        self.embedding[where].insert(0, coalescence)
#        self.embedding[where] += [coalescence]
#        self.sort_embedding()
        lineages = self.maximal_gene_lineages(where)
        assert(all([bool(clade not in lineages) for clade in clades]))
        return(coalescence)
    def genetree(self):
        maximal_gene_lineages = list(self.maximal_gene_lineages())
        assert(len(maximal_gene_lineages)==1)
        genetree = Tree(maximal_gene_lineages.pop())
        if any([clade.name is None for clade in genetree.find_clades()]):
            name_clades(genetree)
        return(genetree)
    def coalesce(self, theta):
        for lineage in self.maximal_gene_lineages():
            lineage.branch_length = 0
        for s_clade in self.tree.find_clades(order='postorder'):
            lineages = self.maximal_gene_lineages(s_clade)
            nr_lineages = len(lineages)
            t = 0
            # Begin coalescence iteration
            while nr_lineages > 1: # While condition (1) is not yet met
                if theta==0: # Deal with theta=0 case separately
                    T = 0
                else:
                    rexp = ExponentialDistribution(float(nr_lineages)*float(nr_lineages-1)/(2*theta))
                    T = rexp.sample()
                if t+T < s_clade.branch_length or self.tree.root==s_clade: # If condition (2) is not yet met
                    # Randomly select two lineages to form the coalescence
                    coalescing_lineages = [list(lineages)[i] for i in np.random.choice(nr_lineages,replace=False,size=2)]
                    # Add the elapsed time to gene tree branch lengths
                    for lineage in lineages:
                        lineage.branch_length += T
                    # Update cumulative time
                    t += T
                    # Coalesce
                    coalescence = self.join(coalescing_lineages, s_clade, branch_length=0)
                    # Update lineages to reflect coalescence
                    lineages = self.maximal_gene_lineages(s_clade) 
                    assert(all([bool(lineage not in lineages) for lineage in coalescing_lineages]))
                    #lineages = (lineages-set(coalescing_lineages)).union(coalescence)
                    nr_lineages = len(lineages)
                else:
                    break
            for lineage in lineages:
                # Add to the gene tree branch lengths the time at the end of the population time span when no coalescence occurred
                if lineage.branch_length is None:
                    lineage.branch_length = 0
                lineage.branch_length += s_clade.branch_length - t
        return(self)

def simulate_coalescent(tree, g_terminal_embedding, theta):
    # theta = N*mu
    # Choose the population (species tree clade: s_clade) which is maximal (with respect to the partial order on the rooted tree) among those which have not been selected already (using depth-first search on the tree with postorder). For this population, get the maximal gene tree clades (gene_lineages) embedded higher in the tree than the given population. These are the lineages that can coalesce in the population (lineages).
    # Coalescence iteration: If the number of lineages (nr_lineages) is greater than 1, draw an exponentially distributed r.var. T with rate equal to nr_lineages*(nr_lineages-1)/(2*theta). If T is less than the time span of the population, randomly choose a pair of gene tree clades (gene_lineages) from lineages, and join these clades (Joining clades is the coalescence operation on the clade embedding object which reduces the number of maximal genetree clades in this embedding by one, see EmbeddedGeneForest.join).
    # Repeat coalescence iteration. Let t_0=0. Starting from iteration i=1, at iteration i record the cumulative time as t_i = t_{i-1}+T. Repeat up to but not including the iteration i such that either: (1) t_i+T is greater than the time span of the population (2) or the number of lineages is 1. If either condition is reached continue for the next maximal population which has not been selected.
    # Assigning branch lengths: Every time a coalescence time a coalescence event occurs for a population (s_clade), add the coalescence time T to the branch length of each lineage in the population. If the cumulative time t exceeds the time span of the population, so that no coalescence occurs in the remaining time, add the remaining time (s_clade.branch_length - t) to each lineage.
    import itertools
    for gene_lineage in itertools.chain.from_iterable(g_terminal_embedding.embedding.values()):
        gene_lineage.branch_length = 0
    gene_lineage_embedding = EmbeddedGeneForest(tree, g_terminal_embedding.embedding)
    gene_lineage_embedding.maximal_gene_lineages()
    for s_clade in tree.find_clades(order='postorder'):
        lineages = gene_lineage_embedding.maximal_gene_lineages(s_clade)
        nr_lineages = len(lineages)
        t = 0
        # Begin coalescence iteration
        while nr_lineages>1: # While condition (1) is not yet met
            if theta==0: # Deal with theta=0 case separately
                T = 0
            else:
                rexp = ExponentialDistribution(float(nr_lineages)*float(nr_lineages-1)/(2*theta))
                T = rexp.sample()
            if t+T < s_clade.branch_length or tree.root==s_clade: # If condition (2) is not yet met
                # Randomly select two lineages to form the coalescence
                coalescing_lineages = [list(lineages)[i] for i in np.random.choice(nr_lineages,replace=False,size=2)]
                # Add the elapsed time to gene tree branch lengths
                for lineage in lineages:
                    lineage.branch_length += T
                # Update cumulative time
                t += T
                # Coalesce
                coalescence = gene_lineage_embedding.join(coalescing_lineages, s_clade, branch_length=0)
                # Update lineages to reflect coalescence
                lineages = gene_lineage_embedding.maximal_gene_lineages(s_clade) 
                assert(all([bool(lineage not in lineages) for lineage in coalescing_lineages]))
                #lineages = (lineages-set(coalescing_lineages)).union(coalescence)
                nr_lineages = len(lineages)
            else:
                break
        for lineage in lineages:
            # Add to the gene tree branch lengths the time at the end of the population time span when no coalescence occurred
            lineage.branch_length += s_clade.branch_length - t
    return(gene_lineage_embedding)
