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

def usage():
    import sys
    cmdref = '{:s} [-k k] [-m m] [--genes=genes] [--indelible_model] [--sims=sims] [--rows=rows] [--cols=cols] [--theta=theta] [--a_max=a_max] [--b_max=b_max] -o <db_file>'.format(sys.argv[0])
    print(cmdref)

def main(argv):
    indelible_model = 'JC'
    indelible_model = 'LAV0.01a'
    theta = 0.01
    mu = 1
    #k = (1,2)
    k = (1,2,3,4,5)
    m = 100
    n = 5
    nr_genes = 10
    nr_sims = 1
    nr_rows = 3
    nr_cols = 3
    a_max = 0.74
    b_max = 0.74
    #a_max = 0.3
    #b_max = 0.2
    kmer_methods = ['CoalescentJCNJ', 'CoalescentJCLS', 'JCNJ','dstarNJ','concatdJCNJ']
    #kmer_methods = ['dstarNJ','concatdJCNJ']
    distance_formulas = ['ARS2015', 'alignment_based']
    multiple_alignment_method = 'clustalo'
    alignment_method = 'stretcher'
    N = theta/mu
    db_file = 'db.sql'

    try:
        opts, args = getopt.getopt(argv,"hk:m:n:o:",["indelible_model=","theta=","genes=","sims=","rows=","cols=","a_max=","b_max="])
    except getopt.GetoptError as err:
        # print usage information and exit:
        print(str(err))
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == "--theta":
            theta = float(arg)
            N = theta
        elif opt == "-k":
            k = tuple([int(k_i) for k_i in re.sub("[()]","",arg).split(",")])
        elif opt == "-m":
            m = int(arg)
        elif opt == "-o":
            db_file = arg
        elif opt == "-n":
            n = int(arg)
        elif opt == "--genes":
            nr_genes = int(arg)
        elif opt == "--method":
            method = arg 
        elif opt == "--distance_formula":
            distance_formula = arg
        elif opt == "--sims":
            nr_sims = int(arg)
        elif opt == "--rows":
            nr_rows = int(arg)
        elif opt == "--cols":
            nr_cols = int(arg)
        elif opt == "--a_max":
            a_max = float(arg)
        elif opt == "--b_max":
            b_max = float(arg)
        elif opt == "--indelible_model":
            indelible_model = arg
        elif opt == "--reconstruct_only":
            reconstruct_only = True
    usage()

    import logging
    sqla_logger = logging.getLogger('sqlalchemy.engine.base.Engine')
    sqla_logger.propagate = False
    sqla_logger.addHandler(logging.FileHandler('/tmp/sqla.log'))

    from sqlalchemy import create_engine
    engine = create_engine('sqlite:///{:s}'.format(os.path.abspath(db_file)), echo=True, convert_unicode=True)
    from sqlalchemy.orm import sessionmaker
    Session = sessionmaker(bind=engine)
    session = Session()
    Base.metadata.create_all(engine)


    import resource
    sim_set = None
    print('Simulating sequence data for a {:d}x{:d} Huelsenbeck diagram with {:d} simulations for each tree, and {:d} gene trees per simulation, using the following parameters: theta={:.2e}, m={:d}, indelible_model={:s}\n'.format(nr_rows,nr_cols,nr_sims,nr_genes,theta,m,indelible_model))
    sim_set = HuelsenbeckSimulationSet(rows=nr_rows, cols=nr_cols, nr_sims=nr_sims, theta=theta, indelible_model=indelible_model, genes=nr_genes, m=m, n=n, a_max=a_max, b_max=b_max)
    for method in ['raxml']:
        tree_estimate_set = HuelsenbeckTreeEstimateSet(simulation_set=sim_set, method=method, distance_formula=None, alignment_method=multiple_alignment_method, k=None)
    session.add(sim_set)
    session.add(tree_estimate_set)
    for row,b in enumerate([(row+1)*(b_max)/nr_rows for row in range(nr_rows)]):
        for col,a in enumerate([(col+1)*(a_max)/nr_cols for col in range(nr_cols)]):
            t_a = abs(-3.0/4.0*log(1-4.0/3.0*a)/mu)
            t_b = abs(-3.0/4.0*log(1-4.0/3.0*b)/mu)
            tree = huelsenbeck_tree(t_a,t_b,5)
            tree_newick = ")".join(tree.format('newick').split(")")[:-1])+")"
            print(tree_newick)

            xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))
            #print(','.join([''.join(split[0])+'|'+''.join(split[1]) for split in xtree.get_splits()])+" (t_a,t_b) = ({:f},{:f}): ".format(t_a,t_b))

            species_set = sorted(tree.get_terminals(),key=lambda species: species.name)
            n = len(species_set)
            species_names = [species.name for species in species_set]
            genes = [GeneLineage(name='s{:d}'.format(i)) for i,_ in enumerate(range(len(species_set)))]
            gene_embedding = dict(zip(species_set,[[gene] for gene in genes]))
            for sim_it in range(nr_sims):
                # Create simulation objects
                sim = Simulation(tree=tree_newick, theta=theta, indelible_model=indelible_model, genes=nr_genes, m=m, n=n)
                session.add(sim)
                huel_sim = HuelsenbeckSimulation(simulation_set=sim_set, simulation=sim, row=row, col=col)
                session.add(huel_sim)
                # Prepare kmer distance matrices to be used to compute averages
                kmer_distance_matrices = dict()
                finite_counts_matrices = dict()
                for distance_formula in ['dstar', 'ARS2015']:
                    kmer_distance_matrices[distance_formula] = dict()
                    finite_counts_matrices[distance_formula] = dict()
                    for k_i in k:
                        kmer_distance_matrices[distance_formula][k_i] = zero_distance_matrix(species_names)
                        finite_counts_matrices[distance_formula][k_i] = zero_distance_matrix(species_names)
                # Prepare concatenated sequence object
                sample_ids = [sample.name for sample in itertools.chain.from_iterable(gene_embedding.values())]
                concatenated_sequences = dict((sample_id,SeqRecord(Seq('',DNAAlphabet()),id=sample_id,name=sample_id,description=sample_id)) for sample_id in sample_ids)

                for gene in range(nr_genes):
                    # generate gene tree and sequences
                    # for each set of sequence, and for each distance formula and value of k, generate a k-mer distance matrix. sum these matrices for all genes
                    # also store the concatenated sequences
                    coalescent = EmbeddedGeneForest(tree, gene_embedding)
                    coalescent.coalesce(theta)
                    genetree = coalescent.genetree()
                    with TemporaryDirectory() as tmpdir:
                        sequences = mutate_indelible(genetree, m, tmpdir, indelible_model, aligned=False)
                        aligned_sequences = SeqIO.to_dict(align_sequences(sequences, multiple_alignment_method, tmpdir))
                        for sample_id in sample_ids:
                            concatenated_sequences[sample_id] += aligned_sequences[sample_id]
                    for distance_formula in ['dstar', 'ARS2015']:
                        for k_i in k:
                            if distance_formula == 'ARS2015':
                                dm = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance, grouping=gene_embedding)
                                #print(dm)
                            elif distance_formula == 'dstar':
                                dm = kmer_distance_matrix(sequences, k_i, dstar_kmer_distance, grouping=gene_embedding)
                            else:
                                raise Exception
                            finite_counts_matrices[distance_formula][k_i] += dm.isfinite()
                            kmer_distance_matrices[distance_formula][k_i] += dm.nantozero()
                            #print(kmer_distance_matrices[distance_formula][k_i])
                for distance_formula in ['dstar', 'ARS2015']:
                    for k_i in k:
                        #print(finite_counts_matrices[distance_formula][k_i])
                        avg_dm = kmer_distance_matrices[distance_formula][k_i]/finite_counts_matrices[distance_formula][k_i]
                        #if distance_formula == 'dstar':
                        #    print(finite_counts_matrices[distance_formula][k_i])
                        #    print(avg_dm)
                        kdm = kmer_distance_matrix_from_dm(avg_dm, sim, distance_formula, None, k_i)
                        session.add(kdm)
                jc_dm = kmer_distance_matrix(concatenated_sequences.values(), 1, aligned_kmer_distance, alignment_fn=stretcher_alignment, grouping=gene_embedding)
                #print(jc_dm)
                kdm = kmer_distance_matrix_from_dm(jc_dm, sim, 'concatdJC', alignment_method, 1)
                session.add(kdm)
                # reconstruct from concatenated sequences using raxml
                with TemporaryDirectory() as tmpdir:
                    t0 = time.clock()
                    xtreehat = RAxML(concatenated_sequences.values(), gene_embedding, tmpdir)
                    t1 = time.clock()

                success = int(xtree.displays(xtreehat))
                print(success)
                tree_estimate = TreeEstimate(simulation=sim, method=tree_estimate_set.method, distance_formula=tree_estimate_set.distance_formula, k=tree_estimate_set.k, splits=','.join([''.join(split[0])+'|'+''.join(split[1]) for split in xtreehat.get_splits()]), success=int(xtree.displays(xtreehat)), dt=t1-t0)
                session.add(tree_estimate)
                #session.commit()
                huel_tree_estimate = HuelsenbeckTreeEstimate(tree_estimate_set=tree_estimate_set, tree_estimate=tree_estimate, huelsenbeck_simulation=huel_sim)
            session.add(huel_tree_estimate)

    session.commit()

    # create tree_estimate sets
    for method in kmer_methods:
        if method in ['CoalescentJCNJ', 'CoalescentJCLS', 'JCNJ']:
            distance_formula = 'ARS2015'
            tree_estimate_set = HuelsenbeckTreeEstimateSet(simulation_set=sim_set, method=method, distance_formula=distance_formula, alignment_method=None, k=",".join([str(k_i) for k_i in k]))
        elif method == 'dstarNJ':
            distance_formula = 'dstar'
            tree_estimate_set = HuelsenbeckTreeEstimateSet(simulation_set=sim_set, method=method, distance_formula=distance_formula, alignment_method=None, k=",".join([str(k_i) for k_i in k]))
        elif method == 'concatdJCNJ':
            distance_formula = 'concatdJC'
            #alignment_method = 'clustalo'
            tree_estimate_set = HuelsenbeckTreeEstimateSet(simulation_set=sim_set, method=method, distance_formula=distance_formula, alignment_method=alignment_method, k="1")
        session.add(tree_estimate_set)
        session.commit()

    # fetch tree_estimate sets that do not require full sequence data
    tree_estimate_sets = session.query(HuelsenbeckTreeEstimateSet).\
                                join(HuelsenbeckTreeEstimateSet.simulation_set).\
                                filter(HuelsenbeckTreeEstimateSet.method.in_(kmer_methods)). \
                                filter(HuelsenbeckTreeEstimateSet.simulation_set==sim_set).all()

    # run tree_estimates
    for tree_estimate_set in tree_estimate_sets:
        method = tree_estimate_set.method
        print(method)
        distance_formula = tree_estimate_set.distance_formula
        #alignment_method = tree_estimate_set.alignment_method
        try:
            k = [int(k_i) for k_i in tree_estimate_set.k.split(",")]
        except AttributeError:
            k = None
        for huel_sim in tree_estimate_set.simulation_set.huelsenbeck_simulations:
            sim = huel_sim.simulation
            treedata = sim.tree
            handle = StringIO(treedata)
            #print(handle.read())
            tree = Phylo.read(handle, "newick")
            xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))
            kmer_distance_matrices = dict((kdm.k,kdm.to_dm()) for kdm in sim.kmer_distance_matrices if kdm.k in k and kdm.distance_formula==distance_formula)

            t0 = time.clock()
            if method == 'CoalescentJCLS':
                xtreehat = TreeMinDistanceFromFiveTaxonCoalescentJCExpectedKmerDistanceParameterizationMap(kmer_distance_matrices)
            elif method == 'CoalescentJCNJ':
                xtreehat = NJArgMinSumOfDistancesFromCoalescentJCExpectedKmerPairDistanceParameterizationMap(kmer_distance_matrices)
            elif method == 'JCNJ':
                #for k,dm in kmer_distance_matrices.items():
                #    print dm
                adjusted_distance_matrices = dict((k,JCKmerDistanceMatrixAdjustment(kmer_distance_matrix,k)) for k,kmer_distance_matrix in kmer_distance_matrices.items()) 
                #for k,dm in adjusted_distance_matrices.items():
                #    print dm
                xtreehat = NJ(adjusted_distance_matrices)
            elif method == 'dstarNJ':
                #for _,dm in kmer_distance_matrices.items():
                #    print dm
                xtreehat = NJ(kmer_distance_matrices)
            elif method == 'concatdJCNJ':
                adjusted_distance_matrices = {1:JCKmerDistanceMatrixAdjustment(kmer_distance_matrices[1],1)}
                xtreehat = NJ(adjusted_distance_matrices)
            else:
                raise(Exception)
            t1 = time.clock()

            success = int(xtree.displays(xtreehat))
            print(','.join([''.join(split[0])+'|'+''.join(split[1]) for split in xtree.get_splits()])+" (t_a,t_b) = ({:f},{:f}): ".format(t_a,t_b)+','.join([''.join(split[0])+'|'+''.join(split[1]) for split in xtreehat.get_splits()])+" ({:d})".format(success))
            #print(k)
            tree_estimate = TreeEstimate(simulation=sim, method=method, distance_formula=distance_formula, k=",".join([str(k_i) for k_i in k]), splits=','.join([''.join(split[0])+'|'+''.join(split[1]) for split in xtreehat.get_splits()]), success=int(xtree.displays(xtreehat)), dt=t1-t0)
            session.add(tree_estimate)
            #session.commit()
            huel_tree_estimate = HuelsenbeckTreeEstimate(tree_estimate_set=tree_estimate_set, tree_estimate=tree_estimate, huelsenbeck_simulation=huel_sim)
            session.add(huel_tree_estimate)
    session.commit()

if __name__ == '__main__':
    #unittest.main()
    import cProfile
    #cProfile.run('main(sys.argv[1:])')
    main(sys.argv[1:])
