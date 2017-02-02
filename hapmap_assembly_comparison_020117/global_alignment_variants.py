### Get variants from global alignment of two sequences
### Version 1.0.0: 02/01/2017
### Author: Felix Francis (felixfrancier@gmail.com)

### Requirements
###

############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()

### Import functions
import datetime
import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO                                 ### to read external alignment files
from Bio.Emboss.Applications import NeedleCommandline   ### for Emboss Needle alignment

############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()



############################################################
### INPUT FILES
############################################################

input_file   =   "aln_Tx303_B73v322kb.txt"
input_file   =   "aln_Tx303_B73v322kb_test.txt"
format = "emboss"
# format = "pair"


seq_1 = "./sequences/TA_2016-10-11T10_46_30_113201_1_25376615_22184_PolyMasked.fasta"
seq_2 = "./sequences/Tx303HPLC_BARDCODE_f_2_Tx303HPLC_BARDCODE_r_2merged_reads_assembly.fasta"

same_strand = 0  # (0: not same strand; 1: same strand)

############################################################
### FUNCTIONS
############################################################
### get rev complement of a sequence
def rev_complement(seq):
	seq = seq.upper()
	basecomplement = {'A':'T', 
					  'C':'G', 
					  'G':'C', 
					  'T':'A', 
					  '-':'-', 
					  'N':'N'}
	letters = list(seq)
	letters = [basecomplement[base] for base in letters]
	complement = (''.join(letters))                                                 #gives the complement of the bases in list letters
	return complement[::-1]     


def get_fasta_seq(input_fasta):
    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta") :
        header, sequence = record.id, str(record.seq)
        return sequence

def pairwise_global_alignment(seq_1, seq_2):
    ### process input sequences
    seq_1 = get_fasta_seq(seq_1)
    seq_2 = get_fasta_seq(seq_2)
    if same_strand == 0:
        seq_1 = rev_complement(seq_1)
    seq_1 = "asis:" + seq_1
    seq_2 = "asis:" + seq_2
    ### needle alignment
    needle_cline = NeedleCommandline(asequence=seq_1, bsequence=seq_2, gapopen=16, gapextend=4, outfile="needle.txt")
    stdout, stderr = needle_cline()

############################################################
### CODE
############################################################


if __name__ =='__main__':
    # pairwise_global_alignment(seq_1, seq_2)

    
    
    
    ### read alignment file
    align = AlignIO.read("needle.txt", "emboss")
    # print(align)

    alignment = AlignIO.read(open("needle.txt"), "emboss")
    # for record in alignment :
        # print(record.seq + " " + record.id)

    print alignment[0].seq



############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))  
print "total time = ", total



