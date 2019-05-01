from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix

from Bio import SeqIO
from Bio import Phylo
import subprocess
import os

from Tree import XTree, enumerate_5_taxon_tree_labelings, make_tree

from scipy.special import binom
from scipy.optimize import minimize

import numpy as np
import itertools
from itertools import combinations
from math import exp,log

#method_parameters = {'CoalescentJCNJ': {'distance_formula':'ARS2015'

def f_xtheta(x,theta,k,m):
    #return(2*(m-k+1)*(1-1.0/4**k*sum([binom(k,i)*3*(3*x)**i/(3+8*theta*i) for i in range(k+1)])))
    return((1-1.0/4**k*sum([binom(k,i)*3*(3*x)**i/(3+8*theta*i) for i in range(k+1)])))
    #return(2*(1-1.0/4**k*sum([binom(k,i)*3*(3*x)**i/(3+8*theta*i) for i in range(k+1)]))) # to be consistent with k-mer distance formula ARS2015

def phi(x,theta,k,l,m):
    return(np.array([f_xtheta(x,theta,k,m),f_xtheta(x,theta,l,m)]))

#def f_ttheta(t,theta,k,m):
    #return(2*(m-k+1)*(1-1.0/4**k*sum([binom(k,i)*(3**(i+1)*exp(-4.0/3*t*i))/(3+8*theta*i) for i in range(k+1)])))
def f_ttheta(t,theta,k):
    try:
        return((1-1.0/4**k*sum([binom(k,i)*(3**(i+1)*exp(-4.0/3*t*i))/(3+8*theta*i) for i in range(k+1)])))
        #return(2*(1-1.0/4**k*sum([binom(k,i)*(3**(i+1)*exp(-4.0/3*t*i))/(3+8*theta*i) for i in range(k+1)]))) # factor of 2 for consistency with k-mer distance formula from ARS2015
    except OverflowError:
        print([t,theta,k])
        return(float('inf'))

def phi_t(t,theta,k1,k2):
    return(np.array([f_ttheta(t,theta,k1),f_ttheta(t,theta,k2)]))

def objfn3_t(x,*args):
    dmk1 = args[0]
    dmk2 = args[1]
    k1 = args[2]
    k2 = args[3]
    #m = args[4]
    names = dmk1.names
    n = len(names)
    #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i+1,n)] for i in range(n)]))
    indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
    return(sum([np.linalg.norm(phi_t(x[idx],x[-1],k1,k2)-np.array([dmk1[names[i],names[j]],dmk2[names[i],names[j]]])) for idx,(i,j) in enumerate(indices)]))



#def f_T(T,k,m):
    #return(2*(m-k+1)*(1-((1+3*exp(-4.0/3*T))/4)**k))
def f_T(T,k):
    return((1-((1+3*exp(-4.0/3*T))/4)**k))
#def psi_T(T,k,l,m):
#    return(np.array([f_T(T,k,m),f_T(T,l,m)]))
def psi_T(T,k1,k2):
    return(np.array([f_T(T,k1),f_T(T,k2)]))

def t_a(a,tree):
    s_terminals = tree.get_terminals()
    n = len(s_terminals)
    t = np.zeros((n,n))
    for i,s_terminal1 in s_terminals:
        for j,s_terminal2 in s_terminals:
            t[i,j] = sum([a[idx] for idx in tree.get_path_index(s_terminal1,s_terminal2)])
    return(t)
def objfn1(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    return(np.linalg.norm(phi(x[0],x[1],k,l,m)-np.array([d_k1,d_k2])))

def objfn1_t(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k1 = args[2]
    k2 = args[3]
    #m = args[4]
    #return(np.linalg.norm(phi_t(x[0],x[1],k,l,m)-np.array([d_k1,d_k2])))
    return(np.linalg.norm(phi_t(x[0],x[1],k1,k2)-np.array([d_k1,d_k2])))

def objfn2(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    theta = args[5]
    return(np.linalg.norm(phi(x[0],theta,k,l,m)-np.array([d_k1,d_k2])))

def objfn2_t(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    theta = args[5]
    return(np.linalg.norm(phi_t(x[0],theta,k,l,m)-np.array([d_k1,d_k2])))

def objfn3(x,*args):
    d_k1 = args[0]
    d_k2 = args[1]
    k = args[2]
    l = args[3]
    m = args[4]
    n_g = d_k1.shape[0]
    return(sum([np.linalg.norm(phi(x[idx],x[-1],k,l,m)-np.array([d_k1[i,j],d_k2[i,j]])) for idx,(i,j) in enumerate(combinations(range(n_g),2))]))


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

def estimate_parameters3(kmer_distance_matrices, mu=None):
    # Here, we use the expected k-mer distance formula which does not use the coalescent model to account for variation in divergence time
    k = kmer_distance_matrices.keys()
    assert(len(k)==2)
    names = kmer_distance_matrices.values()[0].names
    n = len(names)
    # Estimate branch lengths and theta
    #that = np.zeros((n,n))
    bnds = tuple([(0, None)]*(n*(n-1)/2) + [(0, None)])
    opt_result = minimize(objfn3_t, tuple([0.5]*(n*(n-1)/2)+[1]), args=(kmer_distance_matrices[k[0]],kmer_distance_matrices[k[1]],k[0],k[1]), bounds=bnds, method='SLSQP')
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

def reconstruct_tree_raxml(sequences, base_embedding, tmpdir, format="fasta"):
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

def max_d(k):
    return(2*(1.0-1.0/4**k))
def t_dJC(d, k):
    #print("{:d}-mer distance: {:f}".format(k,d))
    eps = 1e-10
    return(-3.0/4.0*log(4.0/3.0*(1.0-min(d,max_d(k)-eps)/2.0)**(1.0/float(k))-1.0/3.0))

def reconstruct_tree_NJ(thatdm):
    # Reconstruct tree
    treehat = DistanceTreeConstructor().nj(thatdm)
    xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
    return(xtreehat)

def JCadjustment(kmer_distance_matrix, k):
    matrix = kmer_distance_matrix.matrix
    thatdm = DistanceMatrix(kmer_distance_matrix.names,[[t_dJC(entry, k) for entry in row] for row in matrix])
    return(thatdm)

def most_common_xtree(xtrees):
    grouped_xtrees = itertools.groupby(xtrees, XTree.get_splits)
    xtree = max([list(g) for k,g in grouped_xtrees],key=len)[0]
    return(xtree)

#def reconstruct_tree_JCNJ(kmer_distance_matrices):
#    candidates = []
#    for k in kmer_distance_matrices.keys():
#        thatdm = JCadjustment(kmer_distance_matrices[k])
#        # Reconstruct tree
#        treehat = DistanceTreeConstructor().nj(thatdm)
#        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
#    #    treehat = DistanceTreeConstructor().nj(kmer_distance_matrices.values()[0])
#        candidates.append(xtreehat)
#    return(most_common_xtree(candidates))

def reconstruct_tree_NJ(kmer_distance_matrices, ignore_coalescent=False):
    candidates = []
    for k in kmer_distance_matrices.keys():
        thatdm = kmer_distance_matrices[k]
        treehat = DistanceTreeConstructor().nj(thatdm)
        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
    #    treehat = DistanceTreeConstructor().nj(kmer_distance_matrices.values()[0])
        candidates.append(xtreehat)
    return(most_common_xtree(candidates))

def reconstruct_tree_CoalescentJCNJ(kmer_distance_matrices):
    candidates = []
    k_pairs = combinations(kmer_distance_matrices.keys(),2)
    for k1,k2 in k_pairs:
        thatdm,thetahat = estimate_parameters3({k1:kmer_distance_matrices[k1],k2:kmer_distance_matrices[k2]}, mu=1.0)
        treehat = DistanceTreeConstructor().nj(thatdm)
        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
        candidates.append(xtreehat)
    return(most_common_xtree(candidates))

def reconstruct_tree_CoalescentJCLS(kmer_distance_matrices):
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
    # Once the rows and columns are ordered in this way, we permute the rows and columns according to a given labeling. When we pass this g_kmer_distances object to the optimization routine, we obtain the minimum distance ||phi_T-D|| where T is 12|345, 123|45 and D=[D_{12},D_{13},D_{14},D_{15},D_{23},D_{24},D_{25},D_{34},D_{35},D_{45}]
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
            opt_result = minimize(objfn6_t, tuple([0.5]*8), args=(kmer_distance_matrices[k_i],labeling,k_i), bounds=bnds, method='SLSQP')
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

# Reconstruction procedure based on existance of 5-taxon tree k-mer distance invariants
def psi_a(a,theta,k):
    t = (a[0]+a[1],a[0]+a[2]+a[5],a[0]+a[3]+a[5]+a[6],a[0]+a[4]+a[5]+a[6],
                    a[1]+a[2]+a[5],a[1]+a[3]+a[5]+a[6],a[1]+a[4]+a[5]+a[6],
                                   a[2]+a[3]+a[6],a[2]+a[4]+a[6],
                                                  a[3]+a[4])
    return(np.array([f_ttheta(t_i,theta,k) for t_i in t]))

def objfn6_t(x,*args):
    kmer_distances_matrix = args[0]
    labeling = args[1]
    k = args[2]
    n = len(kmer_distances_matrix.names)
    indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i+1,n)] for i in range(n)]))
    #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
    kmer_distances = [kmer_distances_matrix[labeling[i],labeling[j]] for i,j in indices]
    return(np.linalg.norm(psi_a(x[:-1],x[-1],k)-np.array(kmer_distances)/2))

def estimate_parameters5(g_kmer_distances, k, m, s_terminals):
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
            opt_result = minimize(objfn6_t, tuple([0.5]*8), args=(g_kmer_distances_k,k_i), bounds=bnds)
        
        if obj < min_obj:
            labeling_min = labeling
            min_obj = obj
    treehat = make_tree([s_terminals[i] for i in labeling])
    return(treehat)
