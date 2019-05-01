import os
import tempfile
from KmerCoalescent.Simulation import TemporaryDirectory
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess

from Bio import pairwise2

from Bio.Emboss.Applications import NeedleCommandline
from Bio.Emboss.Applications import WaterCommandline
from Bio.Emboss.Applications import StretcherCommandline
from Bio import AlignIO

def water_alignment(sequence1,sequence2):
    with TemporaryDirectory() as tmpdir:
        water_fname = os.path.join(tmpdir,"alignment.fas")
        SeqIO.write(sequence1,os.path.join(tmpdir,"seq1.fas"),"fasta")
        SeqIO.write(sequence2,os.path.join(tmpdir,"seq2.fas"),"fasta")
        water_cli = WaterCommandline(asequence=os.path.join(tmpdir,"seq1.fas"), \
                               bsequence=os.path.join(tmpdir,"seq2.fas"), \
                               gapopen=10, \
                               gapextend=0.5, \
                               outfile=water_fname)
        water_cli()
        alignment = AlignIO.read(water_fname, "emboss")
    return(alignment[0],alignment[1])

def stretcher_alignment(sequence1,sequence2):
    with TemporaryDirectory() as tmpdir:
    #tmpdir = tempfile.mkdtemp()
    #print(tmpdir)
        stretcher_fname = os.path.join(tmpdir,"alignment.fas")
        SeqIO.write(sequence1,os.path.join(tmpdir,"seq1.fas"),"fasta")
        SeqIO.write(sequence2,os.path.join(tmpdir,"seq2.fas"),"fasta")
        stretcher_cli = StretcherCommandline(asequence=os.path.join(tmpdir,"seq1.fas"), \
                               bsequence=os.path.join(tmpdir,"seq2.fas"), \
                               gapopen=16, \
                               gapextend=4, \
                               aformat='srspair', \
                               outfile=stretcher_fname)
        stretcher_cli()
        alignment = AlignIO.read(stretcher_fname, "emboss")
    return(alignment[0],alignment[1])

def emboss_alignment(sequence1,sequence2):
    with TemporaryDirectory() as tmpdir:
        needle_fname = os.path.join(tmpdir,"alignment.fas")
        SeqIO.write(sequence1,os.path.join(tmpdir,"seq1.fas"),"fasta")
        SeqIO.write(sequence2,os.path.join(tmpdir,"seq2.fas"),"fasta")
        needle_cli = NeedleCommandline(asequence=os.path.join(tmpdir,"seq1.fas"), \
                               bsequence=os.path.join(tmpdir,"seq2.fas"), \
                               gapopen=10, \
                               gapextend=0.5, \
                               outfile=needle_fname)
        needle_cli()
        alignment = AlignIO.read(needle_fname, "emboss")
    return(alignment[0],alignment[1])

def pairwise2_alignment(sequence1,sequence2):
    alignments = pairwise2.align.globalxx(sequence1, sequence2)
    return(alignments[0][0],alignments[0][1])

def align_sequences(sequences, alignment_method, tmpdir, keep_fasta=False):
    if alignment_method is None:
        return(sequences)
#    with TemporaryDirectory(keep_fasta) as tmpdir:
    if alignment_method=='clustalo':
        nonempty_sequences = list()
        empty_sequences = list()
        for sequence in sequences:
            if len(sequence.seq)>0 and len(sequence.seq)!=sequence.seq.count('-'):
                nonempty_sequences.append(sequence)
            else:
                empty_sequences.append(sequence)
        if len(nonempty_sequences)>1:
            SeqIO.write(nonempty_sequences,os.path.join(tmpdir,"nonempty.fas"),"fasta")
            args = ["clustalo","-i","nonempty.fas","-o","nonempty_alignment.fas"]
            p = subprocess.Popen(args, stdout=subprocess.PIPE, cwd=tmpdir)
            p.wait()
            sequences = list(SeqIO.parse(os.path.join(tmpdir,"nonempty_alignment.fas"), "fasta"))
        elif len(nonempty_sequences)==1:
            sequences = nonempty_sequences
        if len(nonempty_sequences)>0:
            m = len(sequences[0]) # we should have at least one sequence
            for sequence in empty_sequences:
                sequence.seq = Seq("-"*m)
                sequences.append(sequence)
        else: # 0 nonempty sequence
            sequences = empty_sequences
        SeqIO.write(sequences,os.path.join(tmpdir,"alignment.fas"),"fasta")
#    if keep_fasta:
#        return(sequences,os.path.join(tmpdir,"alignment.fas"))
#    else:
    return(sequences)
