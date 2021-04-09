###########################################################################################################
# VPOT - version 2 - 07/01/2021
#                  - added code to handle multiple values in PD and VT fields in TXT file input, when transcript level predictor values are available
#                  - added code for a utility to convert VEP annotated VCF into standard format VCF which can be used by VPOT
# VPOT - version 1 - 26/08/2020
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf, VPOT_1_prioritise, VPOT_2_Gene, VPOT_3_sample_selection, VPOT_4_stats, VPOT_5_merge, VPOT_6_utility #
from shutil import copyfile #
#
#
###########################################################################################################
# Define global variables
##########################################################################################################
suffix=str(int(time.time())) # get a unique timestamp for suffix 
#
supplied_args=0 #
#
tab='\t' # 
nl='\n' #
Maxval=28 # allow for max of 28 characters values for an alpha predictor.
#Maxcoverage=10 # read coverage of sample must be equal or greater then this to included 
Maxcoverage=0 # read coverage of sample must be equal or greater then this to included 
Non_alt_GT_types = ["0","."] #
# for main VPOT prioritisation option 1
info_opt0_msg1=["#VPOT version 2.1 - 09/04/2021 ",
"#tools=$1 # which tool to use -   ",
"#           1: priority - priority tool     ",
"#           2: genef - gene filter       ",
"#           3: samplef - variant filtering ",
"#           4: stats - variant statistics ",
"#           5: merge - merge multiple VPOL output files into one consolidated VPOL file ",
"#           6: utility - VPOT utilities ",
"#inpt1=$2 # for tool ",
"#           1+2+3+4+5 - location for output file+prefix ",
"#              format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/B1",
"#           6 - name of utility function ",
"#              format -  convertVEP -  a utility function to convert VEP annotated VCF into a standard format VCF which can be used by VPOT ",
"#inpt2=$3 # for tool ",
"#           1 - file of input VCF files (1 VCF per line with tab delimiter to sample ID) ",
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/test_inputs/B0_CVM8_split.hg19_multianno.nonintergenic.nonintronic.vcf<tab>SKDP-32.3 ",
"#           2+3+4 - location and name of input post-prioritisation file",
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/final_pV1.txt ",
"#           5 - location of VPOL input file for merge" ,
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/final_pV1.txt ",
"#           6 - depends on utility" ,
"#             - convertVEP - location of VEP annotated VCF" ,
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/VEP_original.vcf ",
"#inpt3=$4 # for tool ",
"#           1 - prioritisation parameters file ",
"#           2 - file of genes for filter/selection - format -  ACTC1 ",
"#           3 - file of samples_IDs with corresponding include/exclude setting - format -  44-1	1 ",
"#                     -  44-1	1  = variant IS in sample 44-1 ",
"#                     -  44-2	0  = variant IS NOT in sample 44-2 ",
"#                     - combination of these values will determine if a variant is maintained.",
"#                     - for above case, a variant is maintain if it is found in 44-1 and not in 44-2. ",
"#                     - Note: if there are more samples than the ones stated, then they do not influence the variant selection. ",
"#           4 - variants in this percentile to include in gene breakdown", #
"#           5 - Not used ", #
"#           6 - depends on utility", #
"#             - convertVEP - name of converted VEP annotated VCF", #
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/VEP_original_converted.vcf ",
"#inpt4=$5 # for tool ",
"#           1 - Not used ",
"#           2 - Not used ",
"#           3 - Sample ID to apply inheritance model filter, if supplied ",
"#           4 - Not used ", #
"#           5 - Not used ", #
"#           6 - depends on utility", #
"#             - convertVEP - not used", #
"#inpt5=$6 # for tool ",
"#           1 - Not used ",
"#           2 - Not used ",
"#           3 - Inheritance model - DN - DeNovo, AD - Autosomal Dominant, AR - Autosomal Recessive, CH - compound Hete, if sample ID supplied ",
"#           4 - Not used ", #
"#           5 - Not used ", #
"#           6 - depends on utility", #
"#             - convertVEP - not used"] #
#
input_type_VCF=True #
sample_loc=-1 #
sample_ID="" #
INFO_loc=-1 #
FORMAT_loc=-1 #
header_ln="" #
blank_variant_ln="" #
#
############################################################################################################
#
###########################################################################################################
def which_option():
	global supplied_args 
	#
#	print "initial_setup():" #
#	print suffix #
#	print sys.argv #
	supplied_args=len(sys.argv) #
#	print supplied_args #
#
	if (supplied_args < 2 ):  # arg [0] is the python program
		for j in range(len(info_opt0_msg1)): #
			print (info_opt0_msg1[j]) #
		return 1 #
	else :
		VPOT_conf.VPOT_option=sys.argv[1] #
#		print (VPOT_conf.VPOT_option) #
	
	return 0 #
#
###########################################################################################################
#
###########################################################################################################
def clean_up():
##
	print ("clean_up():") 
	for f in glob.glob(VPOT_conf.output_dir+"*_"+suffix+"*_tmp.txt") : # 
		os.remove(f) #
#
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	#
	VPOT_conf.init() #
#
	if (which_option() != 0): #
#		print "no good" #
		return #
#
	if (VPOT_conf.VPOT_option=="priority"): #
#		print ("opt1") 
		VPOT_1_prioritise.main() #
	elif (VPOT_conf.VPOT_option=="genef"): #
#		print ("opt2") 
		VPOT_2_Gene.main() #
	elif (VPOT_conf.VPOT_option=="samplef"): #
#		print ("opt3") 
		VPOT_3_sample_selection.main() #
		if (VPOT_conf.inh_model in ["CH"]) : # for CH model we need to only return variants that are for genes that have variants in both parents
#			print ("CH_inh_model") #
			VPOT_conf.input_file=VPOT_conf.temp_output_file #
			VPOT_conf.gene_list=VPOT_conf.sort_file3 #
			VPOT_2_Gene.filter_the_variants() #
	elif (VPOT_conf.VPOT_option=="stats"): #
#		print ("opt4") 
		VPOT_4_stats.main() #
	elif (VPOT_conf.VPOT_option=="merge"): #
#		print ("opt5") 
		VPOT_5_merge.main() #
	elif (VPOT_conf.VPOT_option=="utility"): #
#		print ("opt6") 
		VPOT_6_utility.main() #
	elif (VPOT_conf.VPOT_option=="help"): #
		for j in range(len(info_opt0_msg1)): #
			print (info_opt0_msg1[j]) #
	else :
		print ("VPOT : Invalid tool option specified - please check input command ")
		for j in range(len(info_opt0_msg1)): #
			print (info_opt0_msg1[j]) #
#		
	#clean up working file
	clean_up() #
#
###########################################################################################################
#
# M A I N #
#
# module load python3
# Entry : python3 VPOT_main.py <option> <input file1> <input file2> <input file3>
#
###########################################################################################################
if __name__ == "__main__":  
	print (info_opt0_msg1[0]) #
	main() 
#
