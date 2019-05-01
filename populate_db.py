#%autoindent
#%load_ext autoreload
#%autoreload 2
#%cd /home/cld/Dropbox/multispecies_coalescent/kmers/sim
#%cd /home/cld/Dropbox/multispecies_coalescent/kmers/sim/may16_1753
#%run sim_tree_rec.py
import sys
import os
import tempfile
import shutil
sys.path.append(os.path.join(os.environ["HOME"],"biopython-1.68"))
import unittest
import subprocess

import itertools
from math import log,exp
from copy import copy
from itertools import combinations
from cStringIO import StringIO

import numpy as np
from scipy.special import binom
from scipy.linalg import expm
from scipy.optimize import minimize

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix

from KmerDistanceEvolution.Data import *
import cPickle

indelible_control_template = """
[TYPE] NUCLEOTIDE 1 

[MODEL] JC
  [submodel]    JC
  [indelmodel]   NB  {:f} {:d}  //  Geometric indel length distribution (q=0.4, r=1)
  [insertrate]   {:f}       //  insertion rate = 0.08 relative to substitution rate of 1
  [deleterate]   {:f}       //  deletion rate = 0.12 relative to substitution rate of 1

[MODEL] JCfixed
  [submodel]    JC
  [indelmodel]   NB  0.4 1  //  Geometric indel length distribution (q=0.4, r=1)
  [insertrate]   0.08       //  insertion rate = 0.08 relative to substitution rate of 1
  [deleterate]   0.12       //  deletion rate = 0.12 relative to substitution rate of 1

[TREE] tree {:s};

[PARTITIONS] partition 
  [tree {:s} {:d}]

[EVOLVE] partition {:d} output
"""
#( q, r, # indel model
#  insertionrate,
#  deleterate,
#  tree (e.g. "( (A:1.0, B:1.0):1.0,(C:1.0, D:1.0):1.0)"
#  indelible_model
#  m,
#  N)

def main(argv):
    import numpy as np
    import re
    import sys
#    indelible_model = None
#    theta = 0.01
#    mu = 1
#    #k = (1,2)
#    k = (1,2,3,4,5)
#    m = 100
#    n = 5
#    nr_genes = 100
#    nr_sims = 3
#    nr_rows = 5
#    nr_cols = 5
#    a_max = 0.7
#    b_max = 0.7
#    method = '5taxoninvariant'
#    distance_formulas = ['ARS2015','block']
#    N = theta/mu
#    sim_set_id = None
#    reconstruct_only = False

    import glob
    sim_set_files = glob.glob("sim_set_*.pkl")
    #tree_estimate_set_files = glob.glob("tree_estimate_set_*.pkl")

    import logging
    sqla_logger = logging.getLogger('sqlalchemy.engine.base.Engine')
    sqla_logger.propagate = False
    sqla_logger.addHandler(logging.FileHandler('/tmp/sqla.log'))

    from sqlalchemy import create_engine
    engine = create_engine('sqlite:////home/cld/kmers/db.sql', echo=True)
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    session = Session()
    Base.metadata.create_all(engine)

    import re

    for sim_set_file in sim_set_files:
        m = re.search('sim_set_(.+).pkl', sim_set_file)
        sim_set_idx = m.group(1)
        with open(sim_set_file, "r") as input_file:
            sim_set_dict = cPickle.load(input_file)
        sim_set_attr = sim_set_dict['attr']
        nr_rows = sim_set_attr['rows']
        nr_cols = sim_set_attr['cols']
        nr_sims = sim_set_attr['nr_sims']
        a_max = sim_set_attr['a_max']
        b_max = sim_set_attr['b_max']
        nr_genes = sim_set_attr['genes']
        theta = sim_set_attr['theta']
        m = sim_set_attr['m']
        n = sim_set_attr['n']
        indelible_model = sim_set_attr['indelible_model']
        sim_set = HuelsenbeckSimulationSet(rows=nr_rows, cols=nr_cols, nr_sims=nr_sims, theta=theta, indelible_model=indelible_model, genes=nr_genes, m=m, n=n, a_max=a_max, b_max=b_max)
        #sim_set = HuelsenbeckSimulationSet()
        #for k, v in sim_set_obj['sim_set'].iteritems():
        #        setattr(sim_set, k, v)
        #rows=nr_rows, cols=nr_cols, nr_sims=nr_sims, theta=theta, indelible_model=indelible_model, genes=nr_genes, m=m, n=n, a_max=a_max, b_max=b_max)
        #session.add(sim_set)
        #session.commit()
        ############
        #sim_set = dict(rows=nr_rows, cols=nr_cols, nr_sims=nr_sims, theta=theta, indelible_model=indelible_model, genes=nr_genes, m=m, n=n, a_max=a_max, b_max=b_max)
        ############
        for tree_estimate_set_dict in sim_set_dict['tree_estimate_sets']: 
            method = tree_estimate_set_dict['method']
            alignment_method = tree_estimate_set_dict['alignment_method']
            kvd_method = tree_estimate_set_dict['kvd_method']
            k = tree_estimate_set_dict['k']
            #tree_estimate_set = HuelsenbeckTreeEstimateSet(sim_set_id=sim_set.id, method=method, kvd_method=kvd_method, k=k)
            tree_estimate_set = HuelsenbeckTreeEstimateSet(simulation_set=sim_set, method=method, distance_formula=kvd_method, alignment_method=alignment_method, k=k)
            session.add(tree_estimate_set)
            session.commit()
            tree_estimate_set_dict['id'] = tree_estimate_set.id
        sims = sim_set_dict['sims']
        for sim_container in sims:
            for tree_estimate_container in sim_container['tree_estimates']: 
                sim_dict = sim_container['sim']
                huel_sim_dict = sim_container['huel_sim']
                tree_newick = sim_dict['tree']
                theta = sim_dict['theta']
                indelible_model = sim_dict['indelible_model']
                nr_genes = sim_dict['genes']
                m = sim_dict['m']
                row = huel_sim_dict['row']
                col = huel_sim_dict['col']
                sim = Simulation(tree=tree_newick, theta=theta, indelible_model=indelible_model, genes=nr_genes, m=m)
                session.add(sim)
                #session.commit()
                #huel_sim = HuelsenbeckSimulation(sim_set_id=sim_set.id, sim_id=sim.id, row=row, col=col)
                huel_sim = HuelsenbeckSimulation(simulation_set=sim_set, simulation=sim, row=row, col=col)
                session.add(huel_sim)
                #session.commit()

                tree_estimate_set_id = tree_estimate_container['tree_estimate_set']['id']
                tree_estimate_dict = tree_estimate_container['tree_estimate']
                method = tree_estimate_dict['method']
                kvd_method = tree_estimate_dict['kvd_method']
                #alignment_method = tree_estimate_dict['alignment_method']
                k = tree_estimate_dict['k']
                splits = tree_estimate_dict['splits']
                success = tree_estimate_dict['success']

                #tree_estimate = TreeEstimate(sim_id=sim.id, method=method, kvd_method=kvd_method, k=k, splits=splits, success=success)
                tree_estimate = TreeEstimate(simulation=sim, method=method, distance_formula=kvd_method, k=k, splits=splits, success=success)
                #tree_estimate = TreeEstimate(simulation=sim, method=method, distance_formula=kvd_method, alignment_method=alignment_method, k=k, splits=splits, success=success)
                session.add(tree_estimate)
                #session.commit()
                huel_tree_estimate = HuelsenbeckTreeEstimate(tree_estimate_set_id=tree_estimate_set_id, tree_estimate=tree_estimate, huelsenbeck_simulation=huel_sim)
                session.add(huel_tree_estimate)
    session.commit()
    session.close()

import sys, getopt
if __name__ == '__main__':
    #unittest.main()
    main(sys.argv[1:])
