### adjust the coordinate of a vcf file, so that they are cooresponding to the true coordinates of the reference genome, rather than the alignment.
### Version 1.0.0: 02/13/2017
### Authors: Felix Francis (felixfrancier@gmail.com) 

### Requirements
### vcf files should follow the format described here: https://samtools.github.io/hts-specs/VCFv4.2.pdf



### Import functions
import datetime
import pandas as pd


############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
### INPUTS
############################################################


### line name
reference	= 'CHR1:25375814..25398735'

### input vcf file
input_file		=	"Chr1_25375814_25398735_Tx303HPLC.vcf"					### msa2vcf output file
# input_file		=	"Chr1_25375814_25398735_CML277_HPLC.vcf"					### msa2vcf output file
# input_file		=	"Chr1_25375814_25398735_Mo17_HPLC.vcf"					### msa2vcf output file
# input_file		=	"Chr1_25375814_25398735_Hp301_HPLC.vcf"					### msa2vcf output file
# input_file		=	"Chr1_25375814_25398735_P39_HPLC.vcf"					### msa2vcf output file
input_path		=	"./"
# output_file		=	"coord_adj_Chr1_25375814_25398735_Tx303HPLC.vcf"
output_file		=	"coord_adj_" + input_file

reference_start = 25375814		### resequenced extended coordinates
reference_stop = 25398735		### resequenced extended coordinates


output_file_name = "reference_adj_"+input_file


############################################################
### FUNCTIONS
############################################################

### read first N_lines of vcf file to get the positon of specified lines
def get_line_column_no(selected_line, N_lines):
    with open(input_path +input_file) as myfile:
        lines = [next(myfile) for x in xrange(N_lines)]
        for line in lines:
            line = line.strip("\n").split("\t")
            if line[0] == '#CHROM':
            # if line[0] == 'CHROM':
                headers = line
                column_pos = [i for i, x in enumerate(headers) if x == selected_line][0]
        return int(column_pos)


def adjust_vcf_coordinates(input_path, input_file, output_file, reference):
	coord_adj = 0
	selected_line_column = get_line_column_no(reference, 30)
	with open(input_path + output_file, 'w') as output:	
		with open(input_path +input_file) as f:
			lines = f.readlines()
			for line in lines:
				if line.startswith("#"):
					output.write(line)
				else:
					line = line.strip("\n").split("\t")
					CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, REFERENCE_GENOTYPE, TARGET_GENOTYPE = line
					REFERENCE_GENOTYPE = REFERENCE_GENOTYPE.split(":")[0]
					### make reference genotype 0/0
					if REFERENCE_GENOTYPE == '1/1':
						adj_REF, adj_ALT = ALT, REF
					else:
						adj_REF, adj_ALT = REF, ALT
					# if len(adj_ALT) > len(adj_REF):
						# adj_ALT = '<INS>'
					# elif len(adj_REF) > len(adj_ALT):
						# adj_ALT = '<DEL>'
						
					### adjust coordiate by accounting for reference gaps in alignment for previous genotype pos
					POS = (int(POS) - coord_adj) + reference_start - 1
					
					if len(adj_ALT) > len(adj_REF):
						output.write(str(CHROM)+'\t'+str(POS)+'\t'+str(ID)+'\t'+str(adj_REF)+'\t'+'<INS>'+'\t'+str(QUAL)+'\t'+str(FILTER)+'\t'+str(INFO)+'\t'+str(FORMAT)+'\t'+ '0/0:1' +'\t'+'1/1:1'+ '\n')
					elif len(adj_REF) > len(adj_ALT):
						output.write(str(CHROM)+'\t'+str(POS)+'\t'+str(ID)+'\t'+str(adj_REF)+'\t'+'<DEL>'+'\t'+str(QUAL)+'\t'+str(FILTER)+'\t'+str(INFO)+'\t'+str(FORMAT)+'\t'+ '0/0:1' +'\t'+'1/1:1'+ '\n')
					else:
						output.write(str(CHROM)+'\t'+str(POS)+'\t'+str(ID)+'\t'+str(adj_REF)+'\t'+str(adj_ALT)+'\t'+str(QUAL)+'\t'+str(FILTER)+'\t'+str(INFO)+'\t'+str(FORMAT)+'\t'+ '0/0:1' +'\t'+'1/1:1'+ '\n')
					### adjust coordiate by accounting for reference gaps in alignment for next genotype pos
					if len(adj_ALT) > len(adj_REF):
						# coord_adj += len(adj_ALT) - 1
						coord_adj += len(adj_ALT) - len(adj_REF)
    
############################################################
### CODE
############################################################
        
if __name__ =='__main__':

    ### adjust coordiate by accounting for reference gaps in alignment
	adjust_vcf_coordinates(input_path, input_file, output_file, reference)






    
    
    
    

