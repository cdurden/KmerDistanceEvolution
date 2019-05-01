from ..Util import name_clades
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
    """ This is the precursor to a gene tree, resulting from a (possibly incomplete) sequence of coalescence events. """
    def __init__(self, tree, gene_clade_embedding):
        self.tree = tree
        self.s_terminals = tree.get_terminals()
        assert(all([k in list(tree.find_clades()) for k in gene_clade_embedding.keys()]))
        self.gene_clade_embedding = dict((k,copy(v)) for k,v in gene_clade_embedding.items())
        for species_clade in tree.find_clades(order='postorder'):
            if not self.gene_clade_embedding.has_key(species_clade):
                self.gene_clade_embedding[species_clade] = []
        assert(self.is_sorted())
    def is_sorted(self):
        for species_clade in self.tree.find_clades():
            for i,j in combinations(range(len(self.gene_clade_embedding[species_clade])),2):
                if self.gene_clade_embedding[species_clade][i].is_parent_of(self.gene_clade_embedding[species_clade][j]) and i > j:
                    return False
        return True
    #def maximal_gene_clades(self, population=None):
    def maximal_gene_clades(self, species_clade=None):
        """ Given a population (specified by a node of the species tree), species_clade, (default: root of the tree), return the maximal clades of the gene tree which are embedded at or below that population """
        """ Given a clade of the species tree, clade_s, (default: root of the tree), return the maximal clades of the gene tree which are embedded at or below clade_s """
        # A gene lineage is a node in the gene forest
        # Given an species_clade (default: root of the tree), return the maximal gene_clades which are embedded at or below (i.e. relative to the inverse partial order on the tree, where the root is the maximal element) the given species_clade
        if species_clade is None:
            species_clade = self.tree.root
        maximal_gene_clades = set()
        for species_subclade in species_clade.find_clades(order='preorder'):
            for gene_clade in self.gene_clade_embedding[species_subclade]:
                # We need the lists of clades that comprise the gene_clade_embedding to be sorted so that gene_clades which are maximal, among those embedded within a given species_clade, are reached first. This should be the case because in the join method, we add the clade representing the coalescence to the gene_clade_embedding at index 0, before each of the coalescing clades.
                if not any([bool(maximal_gene_clade.is_parent_of(gene_clade)) for maximal_gene_clade in maximal_gene_clades]):
                    maximal_gene_clades.add(gene_clade)
        return(maximal_gene_clades)
    #def where(self, gene_clade):
    def where(self, gene_clade):
        for i,gene_clades in enumerate(self.gene_clade_embedding.values()):
            if gene_clade in gene_clades:
                return(self.gene_clade_embedding.keys()[i])
    def join(self, clades, where, branch_length=None):
        if branch_length is not None:
            coalescence = Clade(clades=clades, branch_length=branch_length)
        else:
            coalescence = Clade(clades=clades)
        self.gene_clade_embedding[where].insert(0, coalescence)
#        self.gene_clade_embedding[where] += [coalescence]
#        self.sort_gene_clade_embedding()
        maximal_gene_clades = self.maximal_gene_clades(where)
        assert(all([bool(clade not in maximal_gene_clades) for clade in clades]))
        return(coalescence)
    def genetree(self):
        maximal_gene_clades = list(self.maximal_gene_clades())
        assert(len(maximal_gene_clades)==1)
        genetree = Tree(maximal_gene_clades.pop())
        if any([clade.name is None for clade in genetree.find_clades()]):
            name_clades(genetree)
        return(genetree)
    def coalesce(self, theta):
        for maximal_gene_clade in self.maximal_gene_clades():
            maximal_gene_clade.branch_length = 0
        for species_clade in self.tree.find_clades(order='postorder'):
            maximal_gene_clades = self.maximal_gene_clades(species_clade)
            nr_maximal_gene_clades = len(maximal_gene_clades)
            t = 0
            # Begin coalescence iteration
            while nr_maximal_gene_clades > 1: # While condition (1) is not yet met
                if theta==0: # Deal with theta=0 case separately
                    T = 0
                else:
                    rexp = ExponentialDistribution(float(nr_maximal_gene_clades)*float(nr_maximal_gene_clades-1)/(2*theta))
                    T = rexp.sample()
                if t+T < species_clade.branch_length or self.tree.root==species_clade: # If condition (2) is not yet met
                    # Randomly select two maximal_gene_clades to form the coalescence
                    coalescing_maximal_gene_clades = [list(maximal_gene_clades)[i] for i in np.random.choice(nr_maximal_gene_clades,replace=False,size=2)]
                    # Add the elapsed time to gene tree branch lengths
                    for maximal_gene_clade in maximal_gene_clades:
                        maximal_gene_clade.branch_length += T
                    # Update cumulative time
                    t += T
                    # Coalesce
                    coalescence = self.join(coalescing_maximal_gene_clades, species_clade, branch_length=0)
                    # Update maximal_gene_clades to reflect coalescence
                    maximal_gene_clades = self.maximal_gene_clades(species_clade) 
                    assert(all([bool(maximal_gene_clade not in maximal_gene_clades) for maximal_gene_clade in coalescing_maximal_gene_clades]))
                    #maximal_gene_clades = (maximal_gene_clades-set(coalescing_maximal_gene_clades)).union(coalescence)
                    nr_maximal_gene_clades = len(maximal_gene_clades)
                else:
                    break
            for maximal_gene_clade in maximal_gene_clades:
                # Add to the gene tree branch lengths the time at the end of the population time span when no coalescence occurred
                if maximal_gene_clade.branch_length is None:
                    maximal_gene_clade.branch_length = 0
                maximal_gene_clade.branch_length += species_clade.branch_length - t
        return(self)

#def simulate_coalescent(tree, g_terminal_embedding, theta):
#    # theta = N*mu
#    # Choose the population (species tree clade: species_clade) which is maximal (with respect to the partial order on the rooted tree) among those which have not been selected already (using depth-first search on the tree with postorder). For this population, get the maximal gene tree clades (gene_clades) embedded higher in the tree than the given population. These are the maximal_gene_clades that can coalesce in the population (maximal_gene_clades).
#    # Coalescence iteration: If the number of maximal_gene_clades (nr_maximal_gene_clades) is greater than 1, draw an exponentially distributed r.var. T with rate equal to nr_maximal_gene_clades*(nr_maximal_gene_clades-1)/(2*theta). If T is less than the time span of the population, randomly choose a pair of gene tree clades (gene_clades) from maximal_gene_clades, and join these clades (Joining clades is the coalescence operation on the clade embedding object which reduces the number of maximal genetree clades in this embedding by one, see EmbeddedGeneForest.join).
#    # Repeat coalescence iteration. Let t_0=0. Starting from iteration i=1, at iteration i record the cumulative time as t_i = t_{i-1}+T. Repeat up to but not including the iteration i such that either: (1) t_i+T is greater than the time span of the population (2) or the number of maximal_gene_clades is 1. If either condition is reached continue for the next maximal population which has not been selected.
#    # Assigning branch lengths: Every time a coalescence time a coalescence event occurs for a population (species_clade), add the coalescence time T to the branch length of each maximal_gene_clade in the population. If the cumulative time t exceeds the time span of the population, so that no coalescence occurs in the remaining time, add the remaining time (species_clade.branch_length - t) to each maximal_gene_clade.
#    import itertools
#    for gene_clade in itertools.chain.from_iterable(g_terminal_embedding.embedding.values()):
#        gene_clade.branch_length = 0
#    gene_clade_embedding = EmbeddedGeneForest(tree, g_terminal_embedding.embedding)
#    gene_clade_embedding.maximal_gene_clades()
#    for species_clade in tree.find_clades(order='postorder'):
#        maximal_gene_clades = gene_clade_embedding.maximal_gene_clades(species_clade)
#        nr_maximal_gene_clades = len(maximal_gene_clades)
#        t = 0
#        # Begin coalescence iteration
#        while nr_maximal_gene_clades>1: # While condition (1) is not yet met
#            if theta==0: # Deal with theta=0 case separately
#                T = 0
#            else:
#                rexp = ExponentialDistribution(float(nr_maximal_gene_clades)*float(nr_maximal_gene_clades-1)/(2*theta))
#                T = rexp.sample()
#            if t+T < species_clade.branch_length or tree.root==species_clade: # If condition (2) is not yet met
#                # Randomly select two maximal_gene_clades to form the coalescence
#                coalescing_maximal_gene_clades = [list(maximal_gene_clades)[i] for i in np.random.choice(nr_maximal_gene_clades,replace=False,size=2)]
#                # Add the elapsed time to gene tree branch lengths
#                for maximal_gene_clade in maximal_gene_clades:
#                    maximal_gene_clade.branch_length += T
#                # Update cumulative time
#                t += T
#                # Coalesce
#                coalescence = gene_clade_embedding.join(coalescing_maximal_gene_clades, species_clade, branch_length=0)
#                # Update maximal_gene_clades to reflect coalescence
#                maximal_gene_clades = gene_clade_embedding.maximal_gene_clades(species_clade) 
#                assert(all([bool(maximal_gene_clade not in maximal_gene_clades) for maximal_gene_clade in coalescing_maximal_gene_clades]))
#                #maximal_gene_clades = (maximal_gene_clades-set(coalescing_maximal_gene_clades)).union(coalescence)
#                nr_maximal_gene_clades = len(maximal_gene_clades)
#            else:
#                break
#        for maximal_gene_clade in maximal_gene_clades:
#            # Add to the gene tree branch lengths the time at the end of the population time span when no coalescence occurred
#            maximal_gene_clade.branch_length += species_clade.branch_length - t
#    return(gene_clade_embedding)
