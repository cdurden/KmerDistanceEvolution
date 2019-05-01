from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix

from Bio import SeqIO
from Bio import Phylo
import subprocess
import os

from .. import XTree
from ..Util import enumerate_5_taxon_tree_labelings, make_tree
from .DivergenceTime import DistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap, ArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap

from scipy.special import binom
from scipy.optimize import minimize

import numpy as np
import itertools
from itertools import combinations
from math import exp,log

def RAxML(sequences, base_embedding, tmpdir, format="fasta"):
    for sequence in sequences:
        species = [species for species,sample_ids in base_embedding.items() if sequence.id in [sample_id.name for sample_id in sample_ids]][0]
        sequence.id = sequence.name = sequence.description = species.name
    filename = "alignment.fas"
    SeqIO.write(sequences,os.path.join(tmpdir,filename),"fasta")
    #wd,filename = os.path.split(sequence_file_path)
    """ raxmlHPC-AVX -s alignment.fas -m GTRGAMMA --JC69 -p 12345 -n T1 """
    args = ["raxmlHPC-AVX","-s",filename,"-m","GTRGAMMA","--JC69","-p","12345","-n","T1"]
    #p = subprocess.Popen(args, stdout=subprocess.PIPE, cwd=tmpdir)
    #p.wait()
    try:
        subprocess.check_output(args, cwd=tmpdir)
    except subprocess.CalledProcessError as e:
        print(e.output)
    treehat = Phylo.read(os.path.join(tmpdir,"RAxML_bestTree.T1"),"newick")
    xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
    return(xtreehat)

def NJ(thatdm):
    # Reconstruct tree
    treehat = DistanceTreeConstructor().nj(thatdm)
    xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
    return(xtreehat)

def most_common_xtree(xtrees):
    grouped_xtrees = itertools.groupby(xtrees, XTree.get_splits)
    xtree = max([list(g) for k,g in grouped_xtrees],key=len)[0]
    return(xtree)

#def reconstruct_tree_NJ(kmer_distance_matrices, ignore_coalescent=False):
#    candidates = []
#    for k in kmer_distance_matrices.keys():
#        thatdm = kmer_distance_matrices[k]
#        treehat = DistanceTreeConstructor().nj(thatdm)
#        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
#    #    treehat = DistanceTreeConstructor().nj(kmer_distance_matrices.values()[0])
#        candidates.append(xtreehat)
#    return(most_common_xtree(candidates))

def NJArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap(kmer_distance_matrices):
    candidates = []
    k_pairs = combinations(kmer_distance_matrices.keys(),2)
    for k1,k2 in k_pairs:
        thatdm,thetahat = ArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap({k1:kmer_distance_matrices[k1],k2:kmer_distance_matrices[k2]}, mu=1.0)
        treehat = DistanceTreeConstructor().nj(thatdm)
        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
        candidates.append(xtreehat)
    return(most_common_xtree(candidates))

def TreeMinDistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(kmer_distance_matrices):
    # Given g_kmer_distances, and g_terminals (which labels the rows and columns by leaves of the gene tree)
    # We reconstruct a tree with s_terminals (obtained from keys of g_terminal map) at the leaves
    # We do this by considering all possible labelings on a 5-taxon species tree
    # A given labeling of the 5-taxon species tree is specified by a permutation of [0,1,2,3,4].
    # Each entry of this permutation corresponds to an s_terminal, via the following maps:
    # row/column index of g_kmer_distances -> g_terminal -> s_terminal
    # the first map is determined by g_terminals, the second is defined by g_terminal_map
    # Under the hypothesis that entries of s_terminals correspond to the entries of [0,1,2,3,4] with the same indices, and that the entries of g_terminals map to entries of s_terminals in the same way, we can test whether a given labeling best describes the data by considering whether that labeling minimizes the objective function where the g_kmer_distances (D_{ij} = g_kmer_distances[i-1,j-1]) are obtained as [D_{12},D_{13},D_{14},D_{15},D_{23},D_{24},D_{25},D_{34},D_{35},D_{45}] from the parameterization map with the tree given by 12|345, 123|45.
    # If the correspondence between permutations of [0,1,2,3,4], s_terminals, and g_terminals is not this simple one, we reorder the g_kmer_distances so that the latter condition characterizes when a labeling best describes the data.
    # We can think of this in two steps.
    # First, let g_kmer_distances_s denote the g_kmer_distances matrix, with rows and columns reordered so that the first entry of s_terminals corresponds to the first row/column of g_kmer_distances_s, the second ccorresponds to the second, etc.
    # We obtain this by selecting the rows and columns of g_kmer_distances in the order of s_terminals: The first row and column are obtained from g_kmer_distances by selecting the row/column corresponding to the g_terminal which maps to the first s_terminal, and so forth
    # Once the rows and columns are ordered in this way, we permute the rows and columns according to a given labeling. When we pass this g_kmer_distances object to the optimization routine, we obtain the minimum distance ||AlgebraicCoalescentJCExpectedKmerPairDistanceParameterizationMap_T-D|| where T is 12|345, 123|45 and D=[D_{12},D_{13},D_{14},D_{15},D_{23},D_{24},D_{25},D_{34},D_{35},D_{45}]
    names = kmer_distance_matrices.values()[0].names
    assert(len(names)==5)
    k = kmer_distance_matrices.keys()
    obj_mins = dict([(k_i,float("inf")) for k_i in k])
    labeling_mins = dict([(k_i,None) for k_i in k])
    #obj_mins = [float("inf")]*len(k)
    #labeling_mins = [None]*len(k)
    labelings = enumerate_5_taxon_tree_labelings(names)
    #labeling_obj_mins = [dict(zip(labelings,[float("inf")]*len(labelings)))]*len(k)
    labeling_obj_mins = dict([(k_i,dict(zip(labelings,[float("inf")]*len(labelings)))) for k_i in k])
#    print(labeling_obj_mins)
#    for labeling in labelings:
#        M = np.array([[i]*5 for i in range(1,6)])
#        N = np.array([range(1,6)]*5)
#        N[:,labeling][labeling,:]
#        M[:,labeling][labeling,:]

    bnds = tuple([(0, None)]*8)
    for k_i in k:
        for labeling in labelings:
            opt_result = minimize(DistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap, tuple([0.5]*8), args=(kmer_distance_matrices[k_i],labeling,k_i), bounds=bnds, method='SLSQP')
            obj = opt_result.fun
            #print(obj)
            labeling_obj_mins[k_i][labeling] = obj
            if obj < obj_mins[k_i]:
                labeling_mins[k_i] = labeling
                obj_mins[k_i] = obj
#    print(obj_mins)
    #for k_i in k:
    #    labeling_counts[labeling_mins[k_i]] += 1
    candidates = list(set([v for v in labeling_mins.values() if v is not None]))
    if(len(candidates)==0):
        labeling_min = tuple(np.random.permutation(names))
    else:
        #print(candidates)
        obj = dict(zip(candidates,[0]*len(candidates)))
        for k_i in k:
            for labeling in candidates:
#                print(k_i)
#                print(labeling)
                obj[labeling] += labeling_obj_mins[k_i][labeling]
        labeling_min = min(obj,key=obj.get)
    #print(labeling_min)
    treehat = make_tree(labeling_min)
#    xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals() if clade.name in ['A','B','C','D','E']))
    xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
    return(xtreehat)
