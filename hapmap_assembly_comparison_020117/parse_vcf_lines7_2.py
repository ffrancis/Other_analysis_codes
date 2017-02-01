### vcf_conversion: converts a vcf file to ThermoAlign input format
### Version 1.0.0: 06/28/2016
### Authors: Felix Francis (felixfrancier@gmail.com) 

### Requirements
### vcf files should follow the format described here: https://samtools.github.io/hts-specs/VCFv4.2.pdf
### all vcf input files must be named according to the corresponding fasta input files, e.g. chr1.vcf, chr2.vcf, ... 



### Import functions
import datetime
import os
import pandas as pd
import numpy as np

############################################################
### Time to run the code: start timer
############################################################
import time
t0 = time.time()


############################################################
### INPUTS
############################################################

# input_path   =   "/mnt/data27/ffrancis/HapMap3/all_lines/"
input_path   =   "/mnt/data27/ffrancis/HapMap3/HapMap3_01302017/un_imputed_V4coords/"


# select_line     = "Tx303"
selected_line     = '282set_Tx303'

### actual data
# input_file   =   "c1_hmp31_q30.vcf"
input_file   =   input_path + "hmp321_agpv4_chr1.vcf"

### V3 coords
# locus_start = 25375814
# locus_stop = 25398735


### V4 coords qNLB_1_25722269_22589

locus_start = 25722269 
locus_stop = 25744857

### test data
# input_file   =   "c10_first1000.vcf"
# locus_start =  229011
# locus_stop = 229703


############################################################
### FUNCTIONS
############################################################

### read first N_lines of vcf file to get the positon of specified lines
def get_line_column_no(selected_line, N_lines):
    with open(input_file) as myfile:
        lines = [next(myfile) for x in xrange(N_lines)]
        for line in lines:
            line = line.strip("\n").split("\t")
            if line[0] == '#CHROM':
                headers = line
                column_pos = [i for i, x in enumerate(headers) if x == selected_line][0]
        return column_pos

        

############################################################
### CODE
############################################################
        
        
print get_line_column_no(selected_line, 30)


'''
df = pd.read_csv(input_path + input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 3, 4, 720], header=None)
# print df
# df = df[(df.coordinate >= int(start_pos)) & (df.coordinate <= int(stop_pos))]
df = df[(df[1] >= int(locus_start)) & (df[1] <= int(locus_stop))]


df.columns = ['coordinate', 'REF_allele','ALT_allele', "line"]
df.to_csv('parsed_snpsv4020117' + input_file, sep='\t', encoding='utf-8', index=False)
'''
