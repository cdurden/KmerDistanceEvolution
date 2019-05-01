from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix

from Bio import SeqIO
from Bio import Phylo
import subprocess
import os

from .. import XTree
from ..Util import enumerate_5_taxon_tree_labelings, make_tree
from ..Theory import CoalescentJCExpectedKmerPairDistanceParameterizationMap, FiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap

from scipy.special import binom
from scipy.optimize import minimize

import numpy as np
import itertools
from itertools import combinations
from math import exp,log

# Sum, over all pairs of taxa, the distances between CoalescentJCExpectedKmerPairDistanceParameterizationMap and the pairs of (k1-mer and k2-mer) distances
# N.B.: There is one pair of (k1-mer and k2-mer) distances for each pair of taxa. Every pair of taxa contributes to the objective function the distance from this distance pair to the parameterization map.
def SumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap(x,*args):
#def objfn3_t(x,*args):
    dmk1 = args[0]
    dmk2 = args[1]
    k1 = args[2]
    k2 = args[3]
    #m = args[4]
    names = dmk1.names
    n = len(names)
    #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i+1,n)] for i in range(n)]))
    indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
    return(sum([np.linalg.norm(CoalescentJCExpectedKmerPairDistanceParameterizationMap(x[idx],x[-1],k1,k2)-np.array([dmk1[names[i],names[j]],dmk2[names[i],names[j]]])) for idx,(i,j) in enumerate(indices)]))

def t_a(a,tree):
    s_terminals = tree.get_terminals()
    n = len(s_terminals)
    t = np.zeros((n,n))
    for i,s_terminal1 in s_terminals:
        for j,s_terminal2 in s_terminals:
            t[i,j] = sum([a[idx] for idx in tree.get_path_index(s_terminal1,s_terminal2)])
    return(t)

#def objfn1(x,*args):
#    d_k1 = args[0]
#    d_k2 = args[1]
#    k = args[2]
#    l = args[3]
#    m = args[4]
#    return(np.linalg.norm(AlgebraicCoalescentJCExpectedKmerPairDistanceParameterizationMap(x[0],x[1],k,l,m)-np.array([d_k1,d_k2])))

def objfn1_t(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k1 = args[2]
    k2 = args[3]
    #m = args[4]
    #return(np.linalg.norm(AlgebraicCoalescentJCExpectedKmerPairDistanceParameterizationMap_t(x[0],x[1],k,l,m)-np.array([d_k1,d_k2])))
    return(np.linalg.norm(CoalescentJCExpectedKmerPairDistanceParameterizationMap(x[0],x[1],k1,k2)-np.array([d_k1,d_k2])))

def objfn2(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    theta = args[5]
    return(np.linalg.norm(AlgebraicCoalescentJCExpectedKmerPairDistanceParameterizationMap(x[0],theta,k,l,m)-np.array([d_k1,d_k2])))

def objfn2_t(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    theta = args[5]
    return(np.linalg.norm(CoalescentJCExpectedKmerPairDistanceParameterizationMap(x[0],theta,k,l,m)-np.array([d_k1,d_k2])))

def objfn3(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    n_g = d_k1.shape[0]
    return(sum([np.linalg.norm(AlgebraicCoalescentJCExpectedKmerPairDistanceParameterizationMap(x[idx],x[-1],k,l,m)-np.array([d_k1[i,j],d_k2[i,j]])) for idx,(i,j) in enumerate(combinations(range(n_g),2))]))


def objfn4_T(x,*args):
    dmk1 = args[0]
    dmk2 = args[1]
    k1 = args[2]
    k2 = args[3]
    #m = args[4]
    #n_g = d_k1.shape[0]
    #return(sum([np.linalg.norm(psi_T(x[idx],k,l,m)-np.array([d_k[i,j],d_k2[i,j]])) for idx,(i,j) in enumerate(combinations(range(n_g),2))]))
    names = dmk1.names
    n = len(names)
    indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i+1,n)] for i in range(n)]))
    return(sum([np.linalg.norm(psi_T(x[idx],k1,k2)-np.array([dmk1[names[i],names[j]],dmk2[names[i],names[j]]])) for idx,(i,j) in enumerate(indices)]))

def objfn5_T(x,*args):
    d_k1 = args[0]
    k = args[1]
    m = args[2]
    #return((f_T(x,k,m)-d_k)**2)
    return((f_T(x,k)-d_k1)**2)

#def MinimizeSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap(g_kmer_distances, k, m, g_terminal_embedding, mu=None):
def estimate_parameters1(g_kmer_distances, k, m, g_terminal_embedding, mu=None):
    tree = g_terminal_embedding.tree
    s_terminals = tree.get_terminals()
    g_terminals = list(itertools.chain.from_iterable(g_terminal_embedding.embedding[s_terminal] for s_terminal in s_terminals))
    n = len(s_terminals)
    # Estimate branch lengths and theta
    xhat = np.zeros((n,n))
    that = np.zeros((n,n))
    thetahat = np.zeros((n,n))
    #bnds = ((0.00001, 1), (0, None))
    bnds = ((0, None), (0, None))
    for (s_terminal1,s_terminal2) in combinations(s_terminals,2):
#        xhat0 = []
        that0 = []
        thetahat0 = []
        for g_terminal1 in g_terminal_embedding.embedding[s_terminal1]:
            for g_terminal2 in g_terminal_embedding.embedding[s_terminal2]:
                i = g_terminals.index(g_terminal1)
                j = g_terminals.index(g_terminal2)
                opt_result = minimize(objfn1_t, (0,0), args=(g_kmer_distances[0,i,j],g_kmer_distances[1,i,j],k[0],k[1],m), bounds=bnds)
#                xhat0 += [opt_result.x[0]]
                that0 += [opt_result.x[0]]
                thetahat0 += [opt_result.x[1]]
        i = s_terminals.index(s_terminal1)
        j = s_terminals.index(s_terminal2)
        #xhat[i,j] = sum(xhat0)/len(xhat0)
        thetahat[i,j] = sum(thetahat0)/len(thetahat0)
        #that[i,j] = -3.0/4*log(xhat[i,j])
        that[i,j] = sum(that0)/len(that0)
    if mu is not None:
        that /= mu
    return(that, thetahat)

def estimate_parameters2(g_kmer_distances, k, m, g_terminal_embedding, thetahat=None, mu=None):
    tree = g_terminal_embedding.tree
    s_terminals = tree.get_terminals()
    n = len(s_terminals)
    g_terminals = list(itertools.chain.from_iterable(g_terminal_embedding.embedding[s_terminal] for s_terminal in s_terminals))
    if thetahat is None:
        that0, thetahat0 = estimate_parameters1(g_kmer_distances, k, m, g_terminal_embedding)
        thetahat = np.sum(thetahat0)/((n**2-n)/2)
    xhat = np.zeros((n,n))
    that = np.zeros((n,n))
    #bnds = [(0.00001, 1)]
    bnds = [(0, None)]
    for (s_terminal1,s_terminal2) in combinations(s_terminals,2):
        #xhat0 = []
        that0 = []
        thetahat0 = []
        for g_terminal1 in g_terminal_embedding.embedding[s_terminal1]:
            for g_terminal2 in g_terminal_embedding.embedding[s_terminal2]:
                i = g_terminals.index(g_terminal1)
                j = g_terminals.index(g_terminal2)
                opt_result = minimize(objfn2_t, (0,0), args=(g_kmer_distances[0,i,j],g_kmer_distances[1,i,j],k[0],k[1],m), bounds=bnds)
                #xhat0 += [opt_result.x[0]]
                that0 += [opt_result.x[0]]
        i = s_terminals.index(s_terminal1)
        j = s_terminals.index(s_terminal2)
        #xhat[i,j] = sum(xhat0)/len(xhat0)
        #that[i,j] = -3.0/4*log(xhat[i,j])
        that[i,j] = sum(that0)/len(that0)
    if mu is not None:
        that /= mu
    return(that, thetahat)

# 
def ArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap(kmer_distance_matrices, mu=None):
    k = kmer_distance_matrices.keys()
    assert(len(k)==2)
    names = kmer_distance_matrices.values()[0].names
    n = len(names)
    # Estimate branch lengths and theta
    #that = np.zeros((n,n))
    bnds = tuple([(0, None)]*(n*(n-1)/2) + [(0, None)])
    opt_result = minimize(SumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap, tuple([0.5]*(n*(n-1)/2)+[1]), args=(kmer_distance_matrices[k[0]],kmer_distance_matrices[k[1]],k[0],k[1]), bounds=bnds, method='SLSQP')
    thetahat = opt_result.x[-1]
    #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
    that = [[opt_result.x[i*(i-1)/2+j] for j in range(i)]+[0.0] for i in range(n)]
    thatdm = DistanceMatrix(names, that)
#    print dm
#    print thatdm
    return(thatdm,thetahat)

def estimate_parameters4(kmer_distance_matrices, mu=None):
    # Here, we use the expected k-mer distance formula which does not use the coalescent model to account for variation in divergence time
    k = kmer_distance_matrices.keys()
    assert(len(k)==2)
    names = kmer_distance_matrices.values()[0].names
    n = len(names)
    # Estimate branch lengths and theta
    that = np.zeros((n,n))
    bnds = tuple([(0, None)]*((n**2-n)/2))
    opt_result = minimize(objfn4_T, tuple([0.5]*((n**2-n)/2)), args=(kmer_distance_matrices[k[0]],kmer_distance_matrices[k[1]],k[0],k[1]), bounds=bnds, method='SLSQP')
    #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
    that = [[opt_result.x[i*(i-1)/2+j] for j in range(i)]+[0.0] for i in range(n)]
    thatdm = DistanceMatrix(names, that)
    return(thatdm)

def DistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(x,*args):
#def objfn6_t(x,*args):
    kmer_distances_matrix = args[0]
    labeling = args[1]
    k = args[2]
    n = len(kmer_distances_matrix.names)
    indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i+1,n)] for i in range(n)]))
    #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
    kmer_distances = [kmer_distances_matrix[labeling[i],labeling[j]] for i,j in indices]
    return(np.linalg.norm(FiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(x[:-1],x[-1],k)-np.array(kmer_distances)/2))

def MinimizeDistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(g_kmer_distances, k, m, s_terminals):
#def estimate_parameters5(g_kmer_distances, k, m, s_terminals):
    # Here, we use the expected k-mer distance formula which does not use the coalescent model to account for variation in divergence time
    # for each possible 5-taxon tree, minimize ||psi_k,T - d_k||
    # choose the tree which minimizes this norm
    min_obj = float("inf")
    s_terminals = g_terminal_map.keys()
    assert(len(s_terminals)==5)
    S = set([0,1,2,3,4])
    labelings = list()
    for i in S:
        for j in S-set([i,min(S-{i})]):
            labelings.append([i,min(S-{i}),j])
            labelings[-1] += sorted(S-set(labelings[-1]))
    for labeling in labelings:
        M = np.array([[i]*5 for i in range(1,6)])
        N = np.array([range(1,6)]*5)
        N[:,labeling][labeling,:]
        M[:,labeling][labeling,:]

    for labeling in labelings:
        for k_i in k:
            g_kmer_distances_k = g_kmer_distances[k_i][:,labeling][labeling,:]
            opt_result = minimize(DistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap, tuple([0.5]*8), args=(g_kmer_distances_k,k_i), bounds=bnds)
        
        if obj < min_obj:
            labeling_min = labeling
            min_obj = obj
    treehat = make_tree([s_terminals[i] for i in labeling])
    return(treehat)
