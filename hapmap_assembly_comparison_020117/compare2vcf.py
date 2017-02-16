### compare hapmap vcf and resequening vcf
### Version 1.0.0: 02/16/2017
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

### hapmap variants
hapmap_path = "/home/ffrancis/TA_codes_restored/Other_analysis_codes/hapmap_assembly_comparison_020117/HapMap3_perline/"
hapmap_filename = "282set_Tx303_variant_subsetted_merged_flt_c1.vcf"


### resequencing variants
resequencing_path = "/home/ffrancis/TA_codes_restored/Other_analysis_codes/hapmap_assembly_comparison_020117/alignment/"
resequencing_filename = "coord_adj_Chr1_25375814_25398735_Tx303HPLC.vcf"


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
                headers = line
                column_pos = [i for i, x in enumerate(headers) if x == selected_line][0]
        return int(column_pos)



						
############################################################
### CODE
############################################################
        
# if __name__ =='__main__':

    ## adjust coordiate by accounting for reference gaps in alignment
	# adjust_vcf_coordinates(input_path, input_file, output_file, reference)


### hapmap data
hapmap_dict = {}
with open(hapmap_path + hapmap_filename) as hapmap:
	lines = hapmap.readlines()
	for line in lines:
		if line.startswith("Coordinate"):
			pass
		else:
			line = line.strip('\n').split('\t')
			coord, ref, alt, genotype = line
			hapmap_dict[int(coord)] = {'REF':ref.split(','),'ALT':alt.split(','), 'GENOTYPE':genotype}

# print hapmap_dict[25394174]["REF"]
    

### resequence data
resequence_dict = {}
with open(resequencing_path + resequencing_filename) as resequencing:
	lines = resequencing.readlines()
	for line in lines:
		if line.startswith("#"):
			pass
		else:
			line = line.strip('\n').split('\t')
			CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, REFERENCE_GENOTYPE, TARGET_GENOTYPE = line
			resequence_dict[int(POS)] = {'REF':REF.split(','),'ALT':ALT.split(','), 'GENOTYPE':genotype}

# print resequence_dict


### find shared coords
list_coords_hapmap	= list(hapmap_dict.keys())
list_coords_reseq	= list(resequence_dict.keys())
shared_coords = set(list_coords_hapmap).intersection(list_coords_reseq)

no_shared_coords = len(shared_coords)
hapmap_polymorphisms = len(list_coords_hapmap)
reseq_polymorphisms = len(list_coords_reseq)

shared_alleles = 0
shared_zygosity = 0
for coord in shared_coords:
	reseq_ALT = resequence_dict[coord]['ALT']
	reseq_genotype = resequence_dict[coord]['GENOTYPE']
	hapmap_ALT = hapmap_dict[coord]['ALT']
	hapmap_genotype = hapmap_dict[coord]['GENOTYPE']
	
	### select only the homozugous sites
	if len(set(hapmap_genotype.split('/'))) == 1:
		shared_zygosity += 1
		### check if the reseq alt allele is the same as the hapmap alt allele (corresponding to hapmap genotype)
		if reseq_ALT[0] == hapmap_ALT[int(max(hapmap_genotype.split('/')))-1]:
			shared_alleles += 1

	
	

print '# hapmap_polymorphisms = ', hapmap_polymorphisms
print '# reseq_polymorphisms = ', reseq_polymorphisms
print '# no_shared_coords = ', no_shared_coords
print '# shared_zygosity = ', shared_zygosity
print '# shared_alleles = ', shared_alleles	






