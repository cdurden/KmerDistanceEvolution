from Bio.Phylo.TreeConstruction import _DistanceMatrix
import numpy as np
import itertools
from math import log

def zero_distance_matrix(names):
    matrix = [[0 for j in range(i+1)] for i in range(len(names))]
    return(DistanceMatrix(names,matrix))
class DistanceMatrix(_DistanceMatrix):
    def to_np(self):
        names = self.names
        return(np.array([[self[names[i],names[j]] for j in range(len(names))] for i in range(len(names))]))
    def __iadd__(self, x):
        assert(self.names==x.names)
        from operator import add
        for row in range(len(self.matrix)):
            self.matrix[row] = map(add, self.matrix[row], x.matrix[row])
        return(self)
    def __add__(self, x):
        assert(self.names==x.names)
        from operator import add
        matrix = list()
        for row in range(len(self.matrix)):
            matrix.append(map(add, self.matrix[row], x.matrix[row]))
        return(DistanceMatrix(self.names,matrix))
    def __sub__(self, x):
        assert(self.names==x.names)
        from operator import sub
        matrix = list()
        for row in range(len(self.matrix)):
            matrix.append(map(sub, self.matrix[row], x.matrix[row]))
        return(DistanceMatrix(self.names,matrix))
    def __idiv__(self, x):
        assert(self.names==x.names)
        from operator import div
        for row in range(len(self.matrix)):
            self.matrix[row] = map(div, self.matrix[row][:-1], x.matrix[row][:-1]+[0])
        return(self)
    def __div__(self, x):
        assert(self.names==x.names)
        from operator import div 
        matrix = list()
        for row in range(len(self.matrix)):
            row_entries = []
            for col in range(len(self.matrix[row][:-1])):
                try:
                    row_entries.append(self.matrix[row][col]/x.matrix[row][col])
                except ZeroDivisionError:
                    row_entries.append(float('nan'))
            matrix.append(row_entries+[0.0])
            #matrix.append(map(div, self.matrix[row][:-1], x.matrix[row][:-1])+[0])
        return(DistanceMatrix(self.names,matrix))
    def isfinite(self):
        matrix = list()
        for row in range(len(self.matrix)):
            matrix.append(map(lambda x: int(np.isfinite(x)), self.matrix[row]))
        return(DistanceMatrix(self.names,matrix))
    def nantozero(self):
        for row in range(len(self.matrix)):
            self.matrix[row] = [0 if np.isnan(x) else x for x in self.matrix[row]]
        return(self)


def kmer_vector(sequence,k):
    kmer_vector = dict()
    for i in range(len(sequence)-k+1):
        kmer = tuple([sequence[j] for j in range(i,i+k)])
        try: 
            kmer_vector[kmer] = kmer_vector[kmer] + 1
        except KeyError:
            kmer_vector[kmer] = 1
    return(kmer_vector) 
        
#def kmer_distance(kmer_vector_1, kmer_vector_2):
def kmer_distance(sequence1, sequence2, k):
    kv1 = kmer_vector(sequence1,k)
    kv2 = kmer_vector(sequence2,k)
    d = 0 
    for kmer in set(kv1.keys()+kv2.keys()):
        try:
            x = kv1[kmer]
        except KeyError:
            x = 0
        try:
            y = kv2[kmer]
        except KeyError:
            y = 0
        d = d+(x-y)**2
    return(d)

#def scaled_kmer_distance(kmer_vector_1, kmer_vector_2, k, m):
    #m = reduce(lambda x,y: x+y,kmer_vector_1.values())+k-1
    #return(kmer_distance(kmer_vector_1,kmer_vector_2)/float(2*(m-k+1)))
def scaled_kmer_distance(sequence1, sequence2, k):
    assert(len(sequence1)==len(sequence2))
    m = len(sequence1)
    return(kmer_distance(sequence1, sequence2, k)/float(2*(m-k+1)))

def alignment_based_scaled_kmer_distance2(sequence1, sequence2, k):
    nongaps = [sequence1[j]!='-' and sequence2[j]!='-' for j in range(len(sequence1))]
    sequence1 = [sequence1[j] for j,nongap in enumerate(nongaps) if nongap]
    sequence2 = [sequence2[j] for j,nongap in enumerate(nongaps) if nongap]
    m = len(sequence1)
    q = int(floor(m/float(k)))
    d = 0 
    for i in range(q):
        start = k*i
        d += int(any([sequence1[j]!=sequence2[j] for j in range(start,start+k)]))
    d /= float(q)
    return(d)

def aligned_kmer_distance(sequence1, sequence2, k):
    nongaps = [sequence1[j]!='-' and sequence2[j]!='-' for j in range(len(sequence1))]
    sequence1 = [sequence1[j] for j,nongap in enumerate(nongaps) if nongap]
    sequence2 = [sequence2[j] for j,nongap in enumerate(nongaps) if nongap]
    m = len(sequence1)
    if k > m:
        return(float("nan"))
    d = 0
    for i in range(m-k+1):
        d += int(any([sequence1[j]!=sequence2[j] for j in range(i,i+k)]))
    d *= 2.0/float(m-k+1)
    return(d)

def dstar_kmer_distance(sequence1, sequence2, k):
    if '-' in sequence1:
        nongaps = [sequence1[j]!='-' for j in range(len(sequence1))]
        sequence1 = [sequence1[j] for j,nongap in enumerate(nongaps) if nongap]
    if '-' in sequence2:
        nongaps = [sequence2[j]!='-' for j in range(len(sequence2))]
        sequence2 = [sequence2[j] for j,nongap in enumerate(nongaps) if nongap]
    from math import sqrt
    m1 = len(sequence1)
    m2 = len(sequence2)
    if k > min(m1,m2):
        return(float("nan"))
    mu1 = (m1-k+1)/float(4**k)
    mu2 = (m2-k+1)/float(4**k)
    kv1 = kmer_vector(sequence1,k)
    kv2 = kmer_vector(sequence2,k)
    Dstar12 = 0 
    Dstar11 = 0 
    Dstar22 = 0 
#    meanx = 0
    kmers = set(kv1.keys()+kv2.keys())
    #print(m1-k+1)
    #print(m2-k+1)
    try:
        for kmer in kmers:
            try:
                x1 = (kv1[kmer]-mu1)
            except KeyError:
                x1 = -mu1
            try:
                x2 = (kv2[kmer]-mu2)
            except KeyError:
                x2 = -mu2
    #        meanx += x
            Dstar12 += (x1*x2)
            Dstar11 += (x1*x1)
            Dstar22 += (x2*x2)
        Dstar12 += (4**k-len(kmers))*(mu1*mu2)
        Dstar11 += (4**k-len(kmers))*mu1**2
        Dstar22 += (4**k-len(kmers))*mu2**2
        Dstar12 /= sqrt(mu1*mu2)
        Dstar11 /= mu1
        Dstar22 /= mu2
    #    meanx += (4**k-len(kmers))*(-mu1/sqrt(m1-k+1))
        if Dstar12<=0:
            dstar = float('nan')
        else:
            dstar = abs(log(Dstar12/sqrt(Dstar11*Dstar22)))
        return(dstar)
    except ValueError:
        return(float("nan"))
        #print(len(sequence1))
        #print(len(sequence2))
    except ZeroDivisionError:
        return(float("nan"))

def normalized_kmer_distance(sequence1, sequence2, k): #ARS2015
    if '-' in sequence1:
        nongaps = [sequence1[j]!='-' for j in range(len(sequence1))]
        sequence1 = [sequence1[j] for j,nongap in enumerate(nongaps) if nongap]
    if '-' in sequence2:
        nongaps = [sequence2[j]!='-' for j in range(len(sequence2))]
        sequence2 = [sequence2[j] for j,nongap in enumerate(nongaps) if nongap]
    from math import sqrt
    m1 = len(sequence1)
    m2 = len(sequence2)
    if k > min(m1,m2):
        return(float("nan"))
    mu1 = (m1-k+1)/float(4**k)
    mu2 = (m2-k+1)/float(4**k)
    kv1 = kmer_vector(sequence1,k)
    kv2 = kmer_vector(sequence2,k)
    d = 0 
    kmers = set(kv1.keys()+kv2.keys())
    try:
        for kmer in kmers:
            try:
                x = (kv1[kmer]-mu1)/sqrt(m1-k+1)
            except KeyError:
                x = -mu1/sqrt(m1-k+1)
            try:
                y = (kv2[kmer]-mu2)/sqrt(m2-k+1)
            except KeyError:
                y = -mu2/sqrt(m2-k+1)
            d += (x-y)**2
        d += (4**k-len(kmers))*(-mu1/sqrt(m1-k+1)+mu2/sqrt(m2-k+1))**2
        #return(d/2.0) 
        return(d) # This matches ARS2015, need to adjust f_ttheta accordingly
    except ValueError:
        return(float("nan"))
    except ZeroDivisionError:
        return(float("nan"))

def generalized_scaled_kmer_distance(sequence1, sequence2, k):
#    nongaps = [sequence1[j]!='-' and sequence2[j]!='-' for j in range(len(sequence1))]
#    sequence1 = [sequence1[j] for j,nongap in enumerate(nongaps) if nongap]
#    sequence2 = [sequence2[j] for j,nongap in enumerate(nongaps) if nongap]
    m1 = len(sequence1)
    m2 = len(sequence2)
    if m1>m2:
        sequence0 = sequence1
        sequence1 = sequence2
        sequence2 = sequence0
        m0 = m1
        m1 = m2
        m2 = m0
    #kv1 = kmer_vector(sequence1,k)
    #kv2 = kmer_vector([sequence2[i] for i in range(m1)],k)
    #S = scaled_kmer_distance(kv1,kv2,k,m1)
    from math import floor
#    try:
    S = scaled_kmer_distance(sequence1,[sequence2[i] for i in range(m1)],k)
    d = int(floor(m2/float(m1-k+1)))
    for i in range(1,d):
        start = i*(m1-k+1)
        #kv1 = kmer_vector(sequence1,k)
        #kv2 = kmer_vector([sequence2[start+j] for j in range(m1)],k)
        #S += scaled_kmer_distance(kv1,kv2,k,m1)
        S += scaled_kmer_distance(sequence1,[sequence2[start+j] for j in range(m1)],k)
    if m2 % (m1-k+1) > 0:
        start = d*(m1-k+1)-k+1
        #kv1 = kmer_vector([sequence1[i] for i in range(m1-(m2-start),m1)],k)
        #kv2 = kmer_vector([sequence2[i] for i in range(start,m2)],k)
        #T = scaled_kmer_distance(kv1,kv2,k,m2-start)

        #print(start)
        #print(m2)
        #print(m1)
        #print(k)
        T = scaled_kmer_distance([sequence1[i] for i in range(m1-(m2-start),m1)],[sequence2[i] for i in range(start,m2)],k)
        S = (S + (m2-start)/m1*T)/(d+(m2-start)/m1)
    else:
        S = S/d
    return(S)
#    except ZeroDivisionError:
#        return(float("nan"))

def kmer_distance_matrix(sequences, k, kmer_distance_fn, alignment_fn=None, grouping=None):
    n = len(sequences)
    distances = list()
    for i in range(n):
        if alignment_fn is not None:
            row = []
            for j in range(i):
                aligned_sequence1,aligned_sequence2 = alignment_fn(sequences[i],sequences[j])
                row.append(kmer_distance_fn(aligned_sequence1,aligned_sequence2,k))
            distances.append(row+[0.0])
        else:
#            if np.isnan(kmer_distance_fn(sequences[i], sequences[j], k)):
#                print str(sequences[i].seq)
#                print str(sequences[j].seq)
#                print kmer_distance_fn
#                print kmer_distance_fn(sequences[i], sequences[j], k)
            distances.append([kmer_distance_fn(sequences[i], sequences[j], k) for j in range(i)]+[0])
    #print(distances)
    dm = DistanceMatrix([sequence.id for sequence in sequences], matrix=distances)
    #print dm
    if grouping is None:
        return(dm)
    else:
        return(average_over_groups(dm, grouping, sort_by_grouping=True, key=lambda group: group.name))

def average_over_groups(distance_matrix, grouping, sort_by_grouping=True, key=None):
    specieses = grouping.keys()
    if sort_by_grouping:
        specieses.sort(key=key)
    n = len(specieses)
    distances = list()
    for i in range(n):
        distances.append([])
        for j in range(i):
            ungrouped_distances = [distance_matrix[gene1.name,gene2.name] for gene1,gene2 in itertools.product(grouping[specieses[i]],grouping[specieses[j]])]
            distances[i].append(np.nanmean(ungrouped_distances))
        distances[i].append(0)
    return(DistanceMatrix([species.name for species in specieses], matrix=distances))
