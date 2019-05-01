import re
import os, sys, getopt
import itertools
import numpy as np
from math import log
import time

from cStringIO import StringIO
#from StringIO import StringIO
from Bio import Phylo
from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from KmerDistanceEvolution.Data import *
from KmerDistanceEvolution.Data.KmerDistance import normalized_kmer_distance, dstar_kmer_distance, aligned_kmer_distance, kmer_distance_matrix, zero_distance_matrix
from KmerDistanceEvolution import XTree
from KmerDistanceEvolution.Theory import JCKmerDistanceMatrixAdjustment
from KmerDistanceEvolution.Processes.Coalescence import EmbeddedGeneForest
from KmerDistanceEvolution.Processes.Coalescence import GeneLineage
from KmerDistanceEvolution.Processes.Mutation import mutate_indelible
from KmerDistanceEvolution.ParameterEstimation.Tree import TreeMinDistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap, NJArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap, NJ, RAxML
from KmerDistanceEvolution.ParameterEstimation.Alignment import align_sequences, pairwise2_alignment, emboss_alignment, water_alignment, stretcher_alignment
from KmerDistanceEvolution.Util import huelsenbeck_tree, TemporaryDirectory


a = 0.2
b = 0.5
mu = 1
theta = 1
nr_genes = 100
m = 100
k = 2
indelible_model = 'JC'
t_a = abs(-3.0/4.0*log(1-4.0/3.0*a)/mu)
t_b = abs(-3.0/4.0*log(1-4.0/3.0*b)/mu)
tree = huelsenbeck_tree(t_a,t_b,5)
tree_newick = ")".join(tree.format('newick').split(")")[:-1])+")"
print(tree_newick)
xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))

species_set = sorted(tree.get_terminals(),key=lambda species: species.name)
n = len(species_set)
species_names = [species.name for species in species_set]
genes = [GeneLineage(name='s{:d}'.format(i)) for i,_ in enumerate(range(len(species_set)))]
gene_embedding = dict(zip(species_set,[[gene] for gene in genes]))
finite_counts_matrix = zero_distance_matrix(species_names)
kmer_distance_matrix = zero_distance_matrix(species_names)
for gene in range(nr_genes):
    coalescent = EmbeddedGeneForest(tree, gene_embedding)
    coalescent.coalesce(theta)
    genetree = coalescent.genetree()
    with TemporaryDirectory() as tmpdir:
        sequences = mutate_indelible(genetree, m, tmpdir, indelible_model, aligned=False)
    dm = kmer_distance_matrix(sequences, k, normalized_kmer_distance, grouping=gene_embedding)
    finite_counts_matrix += dm.isfinite()
    kmer_distance_matrix += dm.nantozero()

avg_dm = kmer_distance_matrix/finite_counts_matrix
kdm = kmer_distance_matrix_from_dm(avg_dm, sim, distance_formula, None, k_i)
xtreehat_CoalescentJCLS = TreeMinDistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(kmer_distance_matrices)
xtreehat_CoalescentJCNJ = NJArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap(kmer_distance_matrices)
xtreehat_JCNJ = NJ(adjusted_distance_matrices)

if __name__ == '__main__':
    #unittest.main()
    import cProfile
    #cProfile.run('main(sys.argv[1:])')
    main(sys.argv[1:])
