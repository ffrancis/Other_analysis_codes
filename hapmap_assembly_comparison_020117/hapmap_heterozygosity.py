### extracts variants corresponding to a specific line and coordinate range
### Version 1.0.0: 02/01/2017
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

# input_path   =   "/mnt/data27/ffrancis/HapMap3/all_lines/"                                ### first release HapMap3 vcf files V3 coords
input_path   =   "/mnt/data27/ffrancis/HapMap3/HapMap3_01302017/un_imputed_V4coords/"       ### latest HapMap3 vcf files V4 coords

### line name
selected_line     = '282set_Tx303'

### input vcf file
input_file   =    "hmp321_agpv4_chr1.vcf"

### V4 coords qNLB_1_25722269_22589
locus_start = 25722269 
locus_stop = 25744857

### test data
# input_file   =  "test_hmp321_agpv4_chr1.vcf"
# locus_start =  56032
# locus_stop = 57140


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

'''
def subset_hapmap(input_path, input_file,output_file_name):
    selected_line_column = get_line_column_no(selected_line, 30)
    df = pd.read_csv(input_path +input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 3, 4, int(selected_line_column)], header=None)
    df = df[(df[1] >= int(locus_start)) & (df[1] <= int(locus_stop))]
    ### split column by delimiter
    df_column_split = pd.DataFrame(df[selected_line_column].str.split(':',1).tolist(), columns = ['genotype','extra'])
    ### merged the split df with the original df
    df = df.join(df_column_split)
    ### remove redundant columns
    df = df.drop(selected_line_column, 1)
    df = df.drop('extra', 1)
    ### remove columns with no variants called / missing data
    df = df.loc[df['genotype'] != '0/0']
    df = df.loc[df['genotype'] != './.']
    df.columns = ['Coordinate', 'REF_allele','ALT_allele', "Genotype"]
    df.to_csv(output_file_name + input_file[:-4] + ".txt", sep='\t', encoding='utf-8', index=False)
'''

############################################################
### CODE
############################################################
        
# if __name__ =='__main__':
    # subset_hapmap(input_path, input_file,"subsetted_")



selected_line_column = get_line_column_no(selected_line, 30)
df = pd.read_csv(input_path +input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 3, 4, int(selected_line_column)], header=None)
df = df[(df[1] >= int(locus_start)) & (df[1] <= int(locus_stop))]
### split column by delimiter
df_column_split = pd.DataFrame(df[selected_line_column].str.split(':',1).tolist(), columns = ['genotype','extra'])
### merged the split df with the original df
df = df.join(df_column_split)
### remove redundant columns
# df = df.drop(selected_line_column, 1)
df = df.drop('extra', 1)
### remove columns with no variants called / missing data
df = df.loc[df['genotype'] != '0/0']
df = df.loc[df['genotype'] != './.']
df.columns = ['Coordinate', 'REF_allele','ALT_allele', "Selected_line", "Genotype" ]
print df
df.to_csv("subsetted_" + input_file[:-4] + ".txt", sep='\t', encoding='utf-8', index=False)
    
    
    
    
    
    
    
    
    
    

