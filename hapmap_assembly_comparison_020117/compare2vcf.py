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
# hapmap_filename = "282set_Tx303_variant_subsetted_merged_flt_c1.vcf"


### resequencing variants
resequencing_path = "/home/ffrancis/TA_codes_restored/Other_analysis_codes/hapmap_assembly_comparison_020117/alignment/"
# resequencing_filename = "coord_adj_Chr1_25375814_25398735_Tx303_HPLC.vcf"

lines = ['CML277', 'HP301', 'Mo17', 'Tx303', 'P39']

output_file = "output_vcfcomparisons.txt"





############################################################
### FUNCTIONS
############################################################

def compare_vcf_per_line(maize_line):
	hapmap_filename = "282set_" + maize_line + "_variant_subsetted_merged_flt_c1.vcf"
	resequencing_filename = "coord_adj_Chr1_25375814_25398735_" + maize_line + "_HPLC.vcf"

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

	### find shared coords
	list_coords_hapmap	= list(hapmap_dict.keys())
	list_coords_reseq	= list(resequence_dict.keys())
	shared_coords = set(list_coords_hapmap).intersection(list_coords_reseq)

	no_shared_coords = len(shared_coords)
	hapmap_polymorphisms = len(list_coords_hapmap)
	reseq_polymorphisms = len(list_coords_reseq)

	shared_alleles = 0
	shared_zygosity = 0
	shared_coords_final = []
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
				shared_coords_final.append(coord)

	unshared_coords_final = (set(list_coords_hapmap) - set(shared_coords_final))

	unshared_adjacent_indel_counts = 0
	for coord in unshared_coords_final:
		add = False
		if coord - 1 in hapmap_dict:
			if '<INS>' in hapmap_dict[coord - 1]['ALT'] or '<DEL>' in hapmap_dict[coord - 1]['ALT']:
				add = True
		if coord + 1 in hapmap_dict:
			if '<INS>' in hapmap_dict[coord + 1]['ALT'] or '<DEL>' in hapmap_dict[coord + 1]['ALT']:
				add = True
		if add == True:
			unshared_adjacent_indel_counts += 1

	perc_unshared_adj2indels = round(float(unshared_adjacent_indel_counts) / float(len(unshared_coords_final)) *100, 2)
			
	output_data = {'maize_line': maize_line, '#hapmap_polymorphisms' : hapmap_polymorphisms, '#reseq_polymorphisms' : reseq_polymorphisms, 'no_shared_coords' : no_shared_coords, 
	'shared_zygosity' : shared_zygosity, 'shared_alleles' : shared_alleles, '%Unshared_Adj2Indels' : perc_unshared_adj2indels}
	
	return output_data

						
############################################################
### CODE
############################################################
        
if __name__ =='__main__':

    ## compare hapmap and resequencing based variants for each line
	
	with open (output_file, 'w') as output:
		output.write("Line" + "\t" + '#Hapmap_polymorphisms' + "\t" + '#Reseq_polymorphisms' + "\t" + '#Shared_coords' + "\t" 
		+ '#Shared_zygosity' + "\t" + '#Shared_alleles' + "\t" + '%Unshared_Adj2Indels' + "\n" )
		for maize_line in lines:
			output_data = compare_vcf_per_line(maize_line)
			# print output_data
			
			
			Line, Hapmap_polymorphisms, Reseq_polymorphisms, Shared_coords, Shared_zygosity, Shared_alleles, Unshared_Adj2Indels = output_data['maize_line'], output_data['#hapmap_polymorphisms'], output_data['#reseq_polymorphisms'], output_data['no_shared_coords'], output_data['shared_zygosity'], output_data['shared_alleles'], output_data['%Unshared_Adj2Indels']
			
			
			output.write("282_" + Line + "\t" + str(Hapmap_polymorphisms) + "\t" + str(Reseq_polymorphisms) + "\t" + str(Shared_coords) + "\t" + 
			str(Shared_zygosity) + "\t" + str(Shared_alleles) + "\t" + str(Unshared_Adj2Indels) + "\n")





# print output_data

# print '# hapmap_polymorphisms = ', hapmap_polymorphisms
# print '# reseq_polymorphisms = ', reseq_polymorphisms
# print '# no_shared_coords = ', no_shared_coords
# print '# shared_zygosity = ', shared_zygosity
# print '# shared_alleles = ', shared_alleles	
# print '% unshared adj to indels = ', perc_unshared_adj2indels

