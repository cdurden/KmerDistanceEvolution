from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix as DistanceMatrix

from Bio import SeqIO
from Bio import Phylo
import subprocess
import os

from . import XTree
from Util import enumerate_5_taxon_tree_labelings, make_tree

from scipy.special import binom
from scipy.optimize import minimize

import numpy as np
import itertools
from itertools import combinations
from math import exp,log

def AlgebraicCoalescentJCExpectedKmerDistance(x,theta,k):
    #return(2*(m-k+1)*(1-1.0/4**k*sum([binom(k,i)*3*(3*x)**i/(3+8*theta*i) for i in range(k+1)])))
    return((1-1.0/4**k*sum([binom(k,i)*3*(3*x)**i/(3+8*theta*i) for i in range(k+1)])))
    #return(2*(1-1.0/4**k*sum([binom(k,i)*3*(3*x)**i/(3+8*theta*i) for i in range(k+1)]))) # to be consistent with k-mer distance formula ARS2015

def AlgebraicCoalescentJCExpectedKmerPairDistanceParameterizationMap(x,theta,k1,k2):
    return(np.array([AlgebraicCoalescentJCExpectedKmerDistance(x,theta,k1),AlgebraicCoalescentJCExpectedKmerDistance(x,theta,k2)]))

def FiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(a,theta,k):
    t = (a[0]+a[1],a[0]+a[2]+a[5],a[0]+a[3]+a[5]+a[6],a[0]+a[4]+a[5]+a[6],
                    a[1]+a[2]+a[5],a[1]+a[3]+a[5]+a[6],a[1]+a[4]+a[5]+a[6],
                                   a[2]+a[3]+a[6],a[2]+a[4]+a[6],
                                                  a[3]+a[4])
    return(np.array([CoalescentJCExpectedKmerDistance(t_i,theta,k) for t_i in t]))

def CoalescentK80ExpectedKmerDistance(t,theta,kappa,k):
    pass
def CoalescentJCExpectedKmerDistance(t,theta,k):
    try:
        return((1-1.0/4**k*sum([binom(k,i)*(3**(i+1)*exp(-4.0/3*t*i))/(3+8*theta*i) for i in range(k+1)])))
        #return(2*(1-1.0/4**k*sum([binom(k,i)*(3**(i+1)*exp(-4.0/3*t*i))/(3+8*theta*i) for i in range(k+1)]))) # factor of 2 for consistency with k-mer distance formula from ARS2015
    except OverflowError:
        print([t,theta,k])
        return(float('inf'))

def CoalescentJCExpectedKmerPairDistanceParameterizationMap(t,theta,k1,k2):
    return(np.array([CoalescentJCExpectedKmerDistance(t,theta,k1),CoalescentJCExpectedKmerDistance(t,theta,k2)]))

def t_dJC(d, k):
    #print("{:d}-mer distance: {:f}".format(k,d))
    eps = 1e-10
    return(-3.0/4.0*log(4.0/3.0*(1.0-min(d,max_d(k)-eps)/2.0)**(1.0/float(k))-1.0/3.0))

# TreeEstimate procedure based on existance of 5-taxon tree k-mer distance invariants
#def JCadjustment(kmer_distance_matrix, k):
def JCKmerDistanceMatrixAdjustment(kmer_distance_matrix, k):
    matrix = kmer_distance_matrix.matrix
    thatdm = DistanceMatrix(kmer_distance_matrix.names,[[t_dJC(entry, k) for entry in row] for row in matrix])
    return(thatdm)

def f_T(T,k):
    return((1-((1+3*exp(-4.0/3*T))/4)**k))

def psi_T(T,k1,k2):
    return(np.array([f_T(T,k1),f_T(T,k2)]))

def max_d(k):
    return(2*(1.0-1.0/4**k))

