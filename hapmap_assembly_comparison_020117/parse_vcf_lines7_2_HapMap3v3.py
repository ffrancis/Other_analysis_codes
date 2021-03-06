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
# input_path   =   "/mnt/data27/ffrancis/HapMap3/HapMap3_01302017/un_imputed_V4coords/"       ### latest HapMap3 vcf files V4 coords
input_path		=   "/mnt/data27/ffrancis/HapMap3/HapMap3_01302017/un_imputed_V3coords/"       ### latest HapMap3 vcf files V3 coords

### line name
selected_line	= '282set_Tx303'

### input vcf file
# input_file	=    "hmp321_agpv4_chr1.vcf"		### latest HapMap3 vcf files V4 file
input_file		=    "merged_flt_c1.vcf"					### latest HapMap3 vcf files V3 file

### V4 coords qNLB_1_25722269_22589
# locus_start = 25722269 
# locus_stop = 25744857
### V3 coords qNLB_1_25376615_22184
# locus_start = 25376615		### metadata original locus coordinates
# locus_stop = 25398798		### metadata original locus coordinates

locus_start = 25375814		### resequenced extended coordinates
locus_stop = 25398735		### resequenced extended coordinates


### test data
# input_file   =  "test_hmp321_agpv4_chr1.vcf"
# locus_start =  56032
# locus_stop = 57140


output_file_name = "coord_subsettedV3_"


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


def coord_subset_hapmap(input_path, input_file,output_file_name):
    selected_line_column = get_line_column_no(selected_line, 30)
    df = pd.read_csv(input_path +input_file, sep='\t', comment='#', skiprows=0, usecols=[1, 3, 4, int(selected_line_column)], header=None)
    df = df[(df[1] >= int(locus_start)) & (df[1] <= int(locus_stop))]
    df.columns = ['Coordinate', 'REF_allele','ALT_allele', "Genotype"]
    df.to_csv(output_file_name + input_file[:-4] + ".txt", sep='\t', encoding='utf-8', index=False)

def variant_subset_hapmap(output_file_name, input_file):
    df2 = pd.read_csv(output_file_name + input_file[:-4] + ".txt", sep='\t')
    ### split column by delimiter
    df_column_split = pd.DataFrame(df2.Genotype.str.split(':',1).tolist(), columns = ['genotype','extra'])
    ### merged the split df with the original df
    df2 = df2.join(df_column_split)
    ### remove columns with no variants called / missing data
    df2 = df2.loc[df2['genotype'] != '0/0']
    df2 = df2.loc[df2['genotype'] != './.']
    df2 = df2.drop('Genotype', 1)
    df2 = df2.drop('extra', 1)
    df2.to_csv('variant_subsetted_' + input_file[:-4] + '.txt', sep='\t', encoding='utf-8', index=False)
    
############################################################
### CODE
############################################################
        
if __name__ =='__main__':


    ### subset hapmap by coordinates
    coord_subset_hapmap(input_path, input_file,output_file_name)

    
    ### filter based on selected line genotype
    variant_subset_hapmap(output_file_name, input_file)



    
    
    
    

