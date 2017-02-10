### Get variants from multiple sequences
###	a pairwise sequence alignment is carried out, from which the variants are recorded as a vcf file
### Version 1.0.0: 02/09/2017
### Author: Felix Francis (felixfrancier@gmail.com)



############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()

### Import functions
import datetime
import os
import pandas as pd
import subprocess as sp
import numpy as np
import glob
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import AlignIO                                 ### to read external alignment files
# from Bio.Emboss.Applications import NeedleCommandline   ### for Emboss Needle alignment
from Bio.Align.Applications import ClustalOmegaCommandline ### for Emboss ClustalOmega

############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()



############################################################
### INPUT FILES & PATHS
############################################################

msa2vcf_path	= "/home/ffrancis/apps/jvarkit/dist/msa2vcf.jar"
sequence_path	= "../sequences/"
reference_path	= "../sequences/B73_reference/"

format = "emboss"



# reference = "./sequences/TA_2016-10-11T10_46_30_113201_1_25376615_22184_PolyMasked.fasta"		### v3 locus REFERENCE
# reference = "../sequences/v4_locus.fasta"														### v4_locus REFERENCE
reference	= "Chr1_25375814_25398735.fasta"													### v3_locus REFERENCE
query		= "Tx303HPLC.fasta"

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

def pairwise_global_alignment(reference, query):
    ### process input sequences
    reference = get_fasta_seq(reference)
    query = get_fasta_seq(query)
    if same_strand == 0:
        query = rev_complement(query)
    reference = "asis:" + reference
    query = "asis:" + query
    ### needle alignment
    needle_cline = NeedleCommandline(asequence=reference, bsequence=query, gapopen=16, gapextend=4, outfile="needle.txt")
    stdout, stderr = needle_cline()

def clustal_omega_alignment(in_file, out_file):
    ### process input sequences
    # in_file = "asis:" + in_file
    ### needle alignment
    clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True) 
    stdout, stderr = clustalomega_cline()

	
def read_fasta(input_fasta):
	handle	= open(input_fasta, "rU")
	for record in SeqIO.parse(handle, "fasta") :
		header, sequence = record.id, str(record.seq)
		return [header, sequence]
		
############################################################
### CODE
############################################################
if __name__ =='__main__':
	### merge input sequences

	### get all fasta files in sequences diretory
	file_list = glob.glob(str(sequence_path) + "*.fasta")
	for file in file_list:
		query = str(file.split('/')[-1])
		print query
	
		# with open('./merged.fasta', 'w') as outfile:
		with open("./" + str(query[:-6]) + "_merged.fasta", 'w') as outfile:
			r_header_sequence = read_fasta(reference_path + reference)
			r_header, r_sequence = r_header_sequence[0], r_header_sequence[1]
			outfile.write(">"+ r_header + '\n')
			outfile.write(r_sequence + '\n')
			q_header_sequence = read_fasta(sequence_path + query)
			q_header, q_sequence = q_header_sequence[0], q_header_sequence[1]
			q_sequence = rev_complement(q_sequence)
			outfile.write(">"+ q_header + '\n')
			outfile.write(q_sequence + '\n')

		# clustal_omega_alignment("./" + str(query[:-6]) + "_merged.fasta", "clustalo_output")
		clustal_omega_alignment("./" + str(query[:-6]) + "_merged.fasta", str(query[:-6])+ "_alignment_output")


		f0 = open(os.devnull, 'w')
		with open(str(reference[:-6]) + "_" + str(query[:-6]) + ".vcf","wb") as out:
			# sp.Popen(["java", "-jar", "%s" %msa2vcf_path,"./clustalo_output"],stdout=out, stderr=f0)
			sp.Popen(["java", "-jar", "%s" %msa2vcf_path,"./"+str(query[:-6])+ "_alignment_output"],stdout=out, stderr=f0)
		
		
		
############################################################
#Time to run the code: end timer
############################################################
t1 = time.time()
total = t1-t0
total = ("{0:.2f}".format(round(total,2)))  
print "total time = ", total



