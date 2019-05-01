import os
import tempfile
#from . import TemporaryDirectory

indelible_control_template = """
[TYPE] NUCLEOTIDE 1 

[MODEL] model 
{:s}

[TREE] tree {:s};

[PARTITIONS] partition 
  [tree model {:d}]

[EVOLVE] partition {:d} output
"""
#( indelible_model,
#  tree (e.g. "( (A:1.0, B:1.0):1.0,(C:1.0, D:1.0):1.0)"
#  m,
#  N)

indelible_models = {
'JC':
"""
  [submodel]    JC
""",
'LAV0.01a':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.1 100  //  Lavelette (a=1.1, M=100)
  [insertrate]   0.01       // relative to substitution rate of 1
  [deleterate]   0.01
""",
'LAV0.05a':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.1 100 
  [insertrate]   0.05  
  [deleterate]   0.05
""",
'LAV0.1a':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.1 100  
  [insertrate]   0.1       
  [deleterate]   0.1
""",
'LAV0.01b':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.5 100  //  Lavelette (a=1.5, M=100)
  [insertrate]   0.01       // relative to substitution rate of 1
  [deleterate]   0.01
""",
'LAV0.05b':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.5 100 
  [insertrate]   0.05  
  [deleterate]   0.05
""",
'LAV0.1b':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.5 100  
  [insertrate]   0.1       
  [deleterate]   0.1
""",
'LAV0.01c':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.8 100  //  Lavelette (a=1.8, M=100)
  [insertrate]   0.01       // relative to substitution rate of 1
  [deleterate]   0.01
""",
'LAV0.05c':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.8 100 
  [insertrate]   0.05  
  [deleterate]   0.05
""",
'LAV0.1c':
"""
  [submodel]    JC
  [indelmodel]   LAV  1.8 100  
  [insertrate]   0.1       
  [deleterate]   0.1
""",
}

import subprocess

def get_parent(tree, child_clade):
    node_path = tree.get_path(child_clade)
    if len(node_path)==1:
        return(tree.root)
    else:
        return node_path[-2]

        
#def read_alignment():
#    indelible_dir = "/home/cld/pub/research/multispecies_coalescent/indelible/"
#    from Bio import SeqIO
#    records_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(indelible_dir,"alignment.fas"), "fasta"))
#    sequences = SeqIO.parse(os.path.join(indelible_dir,"alignment.fas"), "fasta")

def mutate_indelible(tree, m, tmpdir, indelible_model=None, nodes=None, aligned=False, keep_files=False):
    if nodes is None:
        nodes = tree.get_terminals()
    import re
    tree_newick = re.sub("\)[a-zA-Z0-9]+:","):", ")".join(tree.format('newick').split(")")[:-1])+")")
#    with TemporaryDirectory(keep_files) as tmpdir:
    with open(os.path.join(tmpdir,'control.txt'), 'w') as f:
        q = 0.4
        r = 1
        insertionrate = 0.08
        deleterate = 0.12
        f.write(indelible_control_template.format( indelible_models[indelible_model],
                                                   #")".join(tree.format('newick').split(")")[:-1])+")",
                                                   tree_newick,
                                                   m,
                                                   1))
    args = ["indelible"]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, cwd=tmpdir)
    p.wait()
    from Bio import SeqIO
    from Bio.Seq import Seq
    if aligned:
        #records_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(tmpdir,"output_TRUE.phy"), "phylip-relaxed"))
        sequences = list(SeqIO.parse(os.path.join(tmpdir,"output_TRUE.phy"), "phylip-relaxed"))
    else:
        #records_dict = SeqIO.to_dict(SeqIO.parse(os.path.join(tmpdir,"output.fas"), "fasta"))
        sequences = list(SeqIO.parse(os.path.join(tmpdir,"output.fas"), "fasta"))
#    if keep_files:
#        return(sequences, tmpdir.name)
#    else:
    return(sequences)
