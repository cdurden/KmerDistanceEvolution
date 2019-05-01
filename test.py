import unittest, os, shutil
from Bio import Phylo
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import DNAAlphabet
from cStringIO import StringIO
import numpy as np
from math import log

from KmerDistanceEvolution.Tree import huelsenbeck_tree, XTree
from KmerDistanceEvolution.Simulation.Coalescence import EmbeddedGeneForest
from KmerDistanceEvolution.Simulation.Coalescence import GeneLineage
from KmerDistanceEvolution.Simulation.Mutation import mutate_indelible
from KmerDistanceEvolution.Alignment import align_sequences
from KmerDistanceEvolution.KmerDistance import *
from KmerDistanceEvolution.Estimation import reconstruct_tree_CoalescentJCLS, reconstruct_tree_CoalescentJCNJ, reconstruct_tree_NJ, reconstruct_tree_raxml, JCadjustment
from KmerDistanceEvolution.Estimation import f_ttheta

from KmerDistanceEvolution.Simulation import TemporaryDirectory
from KmerDistanceEvolution.DB import *
class TestCoalescent(unittest.TestCase):
    def test_reconstruct_tree_raxml(self):
        labels = list(np.random.permutation(['A','B','C','D','E']))
        tree = huelsenbeck_tree(-3/4.0*log(0.5),-3/4.0*log(0.8),5)
        k = 2
        theta = 0.01
        m = 100
        specieses = sorted(tree.get_terminals(),key=lambda species: species.name)
        names = [species.name for species in specieses]
        genes = [GeneLineage() for _ in range(len(specieses))]
        base_embedding = dict(zip(specieses,[[gene] for gene in genes]))
        #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
        coalescent = EmbeddedGeneForest(tree, base_embedding)
        coalescent.coalesce(theta)
        genetree = coalescent.genetree()
        #gene_species_map = dict((genes[0].name,species.name) for species,genes in base_embedding.items())
        #concatenated_sequences = dict((species,SeqRecord(Seq('',DNAAlphabet()),id=species,name=species,description=species)) for species in gene_species_map.values())
        sample_ids = [sample.name for sample in itertools.chain.from_iterable(base_embedding.values())]
        concatenated_sequences = dict((sample_id,SeqRecord(Seq('',DNAAlphabet()),id=sample_id,name=sample_id,description=sample_id)) for sample_id in sample_ids)
        with TemporaryDirectory() as tmpdir:
            sequences = mutate_indelible(genetree, m, tmpdir, indelible_model='JC')
            alignment = align_sequences(sequences, 'clustalo', tmpdir, keep_fasta=True)
            for sample_id in sample_ids:
                concatenated_sequences[sample_id] += SeqIO.to_dict(alignment)[sample_id]
            tree = reconstruct_tree_raxml(concatenated_sequences.values(), base_embedding, tmpdir)
#        shutil.rmtree(os.path.dirname(fasta_file))
#        self.assertEqual(str(sequences[zeros[0]].seq),str(sequences[zeros[1]].seq))

    def test_basic_simulation_and_reconstruction_procedure(self):
        from math import log
        # Generate a species tree and a collection of gene lineages embedded at the leaves
        tree = huelsenbeck_tree(-3/4.0*log(0.5),-3/4.0*log(0.8),5)
        specieses = tree.get_terminals()
        genes = [GeneLineage() for _ in range(len(specieses))]
        base_embedding = dict(zip(specieses,[[gene] for gene in genes]))
        # Simulate coalescent process
        coalescent = EmbeddedGeneForest(tree, base_embedding)
        coalescent.coalesce(0.0)
        genetree = coalescent.genetree()
        # Simulate sequences along genetree
        sequences = mutate_indelible(genetree, 10, 'JC', aligned=True)
        aligned_sequences = align_sequences(sequences, 'clustalo')
        # Compute some kmer distances
        k = (1,2,3,4,5)
        k_trunc = (1,3) # for NJ method
        kmer_distance_matrices = dict()
        for k_i in k:
            kmer_distance_matrices[k_i] = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance, grouping=base_embedding)
        # Reconstruct tree 
        treehat1 = reconstruct_tree(dict((k_i,kmer_distance_matrices[k_i]) for k_i in k_trunc))
        treehat2 = reconstruct_tree2(kmer_distance_matrices)
        xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))
        xtreehat1 = XTree(treehat1,dict((clade,set([clade.name])) for clade in treehat1.get_terminals()))
        xtreehat2 = XTree(treehat2,dict((clade,set([clade.name])) for clade in treehat2.get_terminals()))
        xtreehat1.displays(xtree)
        xtreehat2.displays(xtree)

    def test_reconstruct_tree(self):
        # We test the reconstruction procedure by generating distance matrices with entries exactly equal to the expected kmer distances for the given tree.
        # The function f_ttheta maps branch lengths to expected kmer distances.
        labels = list(np.random.permutation(['A','B','C','D','E']))
        a = [np.random.uniform() for i in range(7)]
        #a = [0.6884381370851735, 0.06944367272556096, 0.8481260207499146, 0.02095137803332825, 0.5434220735107143, 0.06683442096191772, 0.7457004343685913] # Tree reconstruction failed for this parameter vector with default optimization method (of scipy.minimize). It succeeds with 'SLSQP' method
        treedata = "((({7:s}:{0:f}, {8:s}:{1:f}):{5:f},{9:s}:{2:f}):0, ({10:s}:{3:f},{11:s}:{4:f}):{6:f}):0".format(*(a+labels))
        handle = StringIO(treedata)
        tree = Phylo.read(handle, "newick")
        k = (1,2)
        theta = 0.01
        specieses = tree.get_terminals()
        names = [species.name for species in specieses]
        #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
        #distances = [[tree.distance(names[i],names[j]) for j in range(i)]+[0.0] for i in range(len(specieses))]
        #dm = DistanceMatrix(names, distances)
        kmer_distance_matrices = dict()
        for k_i in k:
            kmer_distance_matrices[k_i] = DistanceMatrix(names,[[f_ttheta(tree.distance(names[i],names[j]),theta,k_i) for j in range(i)]+[0.0] for i in range(len(specieses))])
        treehat = reconstruct_tree(kmer_distance_matrices)
        xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))
        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
        xtreehat.displays(xtree)
        self.assertTrue(xtreehat.displays(xtree))

    def test_convergence_of_average_kmer_distance(self):
        labels = list(np.random.permutation(['A','B','C','D','E']))
        #a = [np.random.uniform() for i in range(7)]
        #a = [0.6884381370851735, 0.06944367272556096, 0.8481260207499146, 0.02095137803332825, 0.5434220735107143, 0.06683442096191772, 0.7457004343685913] # Tree reconstruction failed for this parameter vector with default optimization method (of scipy.minimize). It succeeds with 'SLSQP' method
        #treedata = "((({7:s}:{0:f}, {8:s}:{1:f}):{5:f},{9:s}:{2:f}):0, ({10:s}:{3:f},{11:s}:{4:f}):{6:f}):0".format(*(a+labels))
        #handle = StringIO(treedata)
        #tree = Phylo.read(handle, "newick")
        tree = huelsenbeck_tree(-3/4.0*log(0.5),-3/4.0*log(0.8),5)
        k = 1
        theta = 0.01
        m = 100
        specieses = sorted(tree.get_terminals(),key=lambda species: species.name)
        names = [species.name for species in specieses]
        genes = [GeneLineage() for _ in range(len(specieses))]
        base_embedding = dict(zip(specieses,[[gene] for gene in genes]))
        #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
        distances = [[tree.distance(names[i],names[j]) for j in range(i)]+[0.0] for i in range(len(specieses))]
        dm = DistanceMatrix(names, distances)
        edm = DistanceMatrix(names,[[f_ttheta(tree.distance(names[i],names[j]),theta,k)*2 for j in range(i)]+[0.0] for i in range(len(specieses))])
        kdm_total = zero_distance_matrix(names)
        finite_counts_matrix = zero_distance_matrix(names)
        N = 1000
        checkpoints = np.linspace(0,N,11)
        for i in range(N):
            coalescent = EmbeddedGeneForest(tree, base_embedding)
            coalescent.coalesce(theta)
            genetree = coalescent.genetree()
            with TemporaryDirectory() as tmpdir:
                sequences = mutate_indelible(genetree, m, tmpdir, indelible_model='JC')
#            self.assertEqual(str(sequences[zeros[0]].seq),str(sequences[zeros[1]].seq))
            #kdm = kmer_distance_matrix(sequences, k, scaled_kmer_distance, grouping=base_embedding)
            kdm = kmer_distance_matrix(sequences, k, normalized_kmer_distance, grouping=base_embedding)
            kdm_total += kdm
            finite_counts_matrix += dm.isfinite()
        adm = kdm_total/finite_counts_matrix
        print edm
        print adm
        print edm-adm
#        print dm
#        print edm
#        print adm
#        indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(len(specieses))]))
#        indices.sort(key=lambda (i,j): dm[i,j])
#        indices
#        indices.sort(key=lambda (i,j): edm[i,j])
#        indices
#        indices.sort(key=lambda (i,j): adm[i,j])
#        indices

 
    def test_reconstruct_tree2(self):
        labels = list(np.random.permutation(['A','B','C','D','E']))
        a = [np.random.uniform() for i in range(7)]
        treedata = "((({7:s}:{0:f}, {8:s}:{1:f}):{5:f},{9:s}:{2:f}):0, ({10:s}:{3:f},{11:s}:{4:f}):{6:f}):0".format(*(a+labels))
        handle = StringIO(treedata)
        tree = Phylo.read(handle, "newick")
        k = (1,2,3,4,5)
        theta = 0.01
        specieses = tree.get_terminals()
        names = [species.name for species in specieses]
        #indices = list(itertools.chain.from_iterable([[(i,j) for j in range(i)] for i in range(n)]))
        #distances = [[tree.distance(names[i],names[j]) for j in range(i)]+[0.0] for i in range(n)]
        #dm = DistanceMatrix(names, distances)
        kmer_distance_matrices = dict()
        for k_i in k:
            kmer_distance_matrices[k_i] = DistanceMatrix(names,[[f_ttheta(tree.distance(names[i],names[j]),theta,k_i) for j in range(i)]+[0.0] for i in range(len(specieses))])
            #kmer_distance_matrices[k_i] = DistanceMatrix(names,[[f_ttheta(tree.distance(names[i],names[j]),theta,k_i) for j in np.random.permutation(range(i))]+[0.0] for i in range(len(specieses))]) # This should fail more often than not
        treehat = reconstruct_tree2(kmer_distance_matrices)
        xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))
        xtreehat = XTree(treehat,dict((clade,set([clade.name])) for clade in treehat.get_terminals()))
        self.assertTrue(xtreehat.displays(xtree))
 
    def test_identity_of_nondivergent_sequences(self):
        labels = ['A','B','C','D','E']
        m = 100
        k = (1,2,3,4,5)
        for zeros in [(0,1),(0,2,5),(1,2,5),(2,3,6),(2,4,6),(3,4)]:
            a = [1 if i not in zeros else 0 for i in range(7)]
            treedata = "((({7:s}:{0:f}, {8:s}:{1:f}):{5:f},{9:s}:{2:f}):0, ({10:s}:{3:f},{11:s}:{4:f}):{6:f}):0".format(*(a+labels))
            handle = StringIO(treedata)
            tree = Phylo.read(handle, "newick")
            sequences = mutate_indelible(tree, m, indelible_model='JC')
            self.assertEqual(str(sequences[zeros[0]].seq),str(sequences[zeros[1]].seq))
            for k_i in k:
                dm = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance)
                self.assertEqual(dm[labels[zeros[0]],labels[zeros[1]]],0)

    def test_identification_of_nondivergent_species(self):
        labels = ['A','B','C','D','E']
        m = 100
        k = (1,2,3,4,5)
        for zeros in [(0,1),(0,2,5),(1,2,5),(2,3,6),(2,4,6),(3,4)]:
            a = [1 if i not in zeros else 0 for i in range(7)]
            treedata = "((({7:s}:{0:f}, {8:s}:{1:f}):{5:f},{9:s}:{2:f}):0, ({10:s}:{3:f},{11:s}:{4:f}):{6:f}):0".format(*(a+labels))
            handle = StringIO(treedata)
            tree = Phylo.read(handle, "newick")
            specieses = sorted(tree.get_terminals(),key=lambda species: species.name)
            genes = [GeneLineage() for _ in range(len(specieses))]
            base_embedding = dict(zip(specieses,[[gene] for gene in genes]))
            coalescent = EmbeddedGeneForest(tree, base_embedding)
            coalescent.coalesce(0)
            genetree = coalescent.genetree()
            sequences = mutate_indelible(genetree, m, indelible_model='JC')
#            self.assertEqual(str(sequences[zeros[0]].seq),str(sequences[zeros[1]].seq))
            for k_i in k:
                dm = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance, grouping=base_embedding)
                self.assertEqual(dm[labels[zeros[0]],labels[zeros[1]]],0)

    def test_kmer_distance_storage_and_retrieval(self):
        from sqlalchemy import create_engine
        with TemporaryDirectory() as tmpdir:
            import logging
            sqla_logger = logging.getLogger('sqlalchemy.engine.base.Engine')
            sqla_logger.propagate = False
            sqla_logger.addHandler(logging.FileHandler('/tmp/sqla.log'))

            engine = create_engine('sqlite:///{:s}'.format(os.path.join(tmpdir,'db.sql')), echo=True, convert_unicode=True)
            from sqlalchemy.orm import sessionmaker
            Session = sessionmaker(bind=engine)
            session = Session()
            Base.metadata.create_all(engine)
    
            labels = ['A','B','C','D','E']
            m = 100
            k = (1,2,3,4,5)
            indelible_model = 'JC'
            nr_genes = 100
            kmer_formula = 'ARS2015'
            alignment_method = None
            a = [np.random.uniform() for i in range(7)]
            treedata = "((({7:s}:{0:f}, {8:s}:{1:f}):{5:f},{9:s}:{2:f}):0, ({10:s}:{3:f},{11:s}:{4:f}):{6:f}):0".format(*(a+labels))
            handle = StringIO(treedata)
            tree = Phylo.read(handle, "newick")
            specieses = sorted(tree.get_terminals(),key=lambda species: species.name)
            genes = [GeneLineage() for _ in range(len(specieses))]
            base_embedding = dict(zip(specieses,[[gene] for gene in genes]))
            coalescent = EmbeddedGeneForest(tree, base_embedding)
            coalescent.coalesce(0)
            genetree = coalescent.genetree()
            sequences = mutate_indelible(genetree, m, indelible_model=indelible_model)
#            self.assertEqual(str(sequences[zeros[0]].seq),str(sequences[zeros[1]].seq))
            tree_newick = ")".join(tree.format('newick').split(")")[:-1])+")"
            sim = Simulation(tree=tree_newick, theta=0, indelible_model=indelible_model, genes=nr_genes, m=m)
            session.add(sim)
            for k_i in k:
                dm = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance, grouping=base_embedding)
                kdm = kmer_distance_matrix_from_dm(dm, sim, kmer_formula, alignment_method, k_i)
                session.add(kdm)
                session.commit()
                kdm2 = [kdm for kdm in sim.kmer_distance_matrices if kdm.k==k_i][0]
#                kdm2 = session.query(KmerDistanceMatrix).\
#                        join(KmerDistanceMatrix.simulation).\
#                        filter(KmerDistanceMatrix.simulation==sim).\
#                        filter(KmerDistanceMatrix.k==k_i).\
#                        filter(KmerDistanceMatrix.kmer_formula==kmer_formula).\
#                        filter(KmerDistanceMatrix.alignment_method==alignment_method).one()
                self.assertEqual(dm.matrix,kdm2.to_dm().matrix)

if __name__ == '__main__':
    unittest.main()

#%cd /home/cld/kmers/sim
#%autoindent
#from Tree import huelsenbeck_tree, XTree
#from Simulation.Coalescent import EmbeddedGeneForest
#from Simulation.Coalescent import GeneLineage
#from Simulation.Mutation import mutate_indelible
#from Alignment import align_sequences
#from KmerDistance import *
#from Estimation import reconstruct_tree, reconstruct_tree2
#
#from math import log
#tree = huelsenbeck_tree(-3/4.0*log(0.5),-3/4.0*log(0.8),4)
#tree = huelsenbeck_tree(-3/4.0*log(0.5),0,5)
#tree = huelsenbeck_tree(0.383119,0.107326,5)
#specieses = tree.get_terminals()
#genes = [GeneLineage() for _ in range(len(specieses))]
#base_embedding = dict(zip(specieses,[[gene] for gene in genes]))
#
#coalescent = EmbeddedGeneForest(tree, base_embedding)
#coalescent.coalesce(0.0)
#genetree = coalescent.genetree()
#print genetree
#
#sequences = mutate_indelible(genetree, 10, 'JC', aligned=True)
#aligned_sequences = align_sequences(sequences, 'clustalo')
#k = (1,2,3,4,5)
#kmer_distance_matrices = dict()
#for k_i in k:
#    #kmer_distance_matrices[k_i] = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance, grouping=base_embedding)
#    kmer_distance_matrices[k_i] = kmer_distance_matrix(sequences, k_i, normalized_kmer_distance)
##    print(kmer_distance_matrices[k_i].to_np())
#
#from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
#treehat3 = DistanceTreeConstructor().nj(kmer_distance_matrices[k_i])
#
#treehat1 = reconstruct_tree(dict((k_i,kmer_distance_matrices[k_i]) for k_i in [1,3]))
#treehat2 = reconstruct_tree2(kmer_distance_matrices)
#
#xtree = XTree(tree,dict((clade,set([clade.name])) for clade in tree.get_terminals()))
#xtreehat1 = XTree(treehat1,dict((clade,set([clade.name])) for clade in treehat1.get_terminals()))
#xtreehat2 = XTree(treehat2,dict((clade,set([clade.name])) for clade in treehat2.get_terminals()))
#xtreehat3 = XTree(treehat3,dict((clade,set([clade.name])) for clade in treehat3.get_terminals()))
#
#xtree.get_splits()
#xtreehat1.get_splits()
#xtreehat2.get_splits()
#xtreehat3.get_splits()
#
#xtreehat1.displays(xtree)
#xtreehat2.displays(xtree)
#xtreehat3.displays(xtree)
#
#from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
#sequences = [
#SeqRecord('AAAAAA','A','A'),
#SeqRecord('TTTTTT','B','B'),
#SeqRecord('AAAAAA','C','C'),
#SeqRecord('TTTAAA','D','D'),
#SeqRecord('TTTCCC','E','E'),
#]
#dm1 = kmer_distance_matrix(sequences,2,normalized_kmer_distance)
#dm2 = kmer_distance_matrix(sequences,3,normalized_kmer_distance)
#dm1.to_np()
#dm2.to_np()
#kmer_distance_matrices = {2: dm1, 3: dm2}
#treehat = reconstruct_tree(kmer_distance_matrices)
#treehat = reconstruct_tree2(kmer_distance_matrices)
#print treehat1
#print treehat2
#print genetree

