###########################################################################################################
# 
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf, VPOT_1_prioritise, VPOT_2_Gene, VPOT_3_sample_selection, VPOT_4_stats, VPOT_5_merge #
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
info_opt0_msg1=["#tools=$1 # which tool to use -   ",
"#           1: priority - priority tool     ",
"#           2: genef - gene filter       ",
"#           3: samplef - variant filtering ",
"#           4: stats - variant statistics ",
"#           5: merge - merge two output files into one consolidated output file ",
"#inpt1=$2 # for tool ",
"#           1+2+3+4+5 - location for output file+prefix ",
"#              format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/B1",
"#inpt2=$3 # for tool ",
"#           1 - file of input VCF files (1 VCF per line with tab delimiter to sample ID) ",
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/test_inputs/B0_CVM8_split.hg19_multianno.nonintergenic.nonintronic.vcf<tab>SKDP-32.3 ",
"#           2+3+4 - location and name of input post-prioritisation file",
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/final_pV1.txt ",
"#           5 - location of one input file for merge" ,
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/final_pV1.txt ",
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
"#           5 - location of second input file for merge" ,
"#               format -  /short/a32/exi569/WGS_model/variant_prioritisation_tool/output/final_pV2.txt "] #
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
	elif (VPOT_conf.VPOT_option=="stats"): #
#		print ("opt4") 
		VPOT_4_stats.main() #
	elif (VPOT_conf.VPOT_option=="merge"): #
#		print ("opt5") 
		VPOT_5_merge.main() #
	else :
		print ("VPOT : Invalid tool option specified - please check input command ")
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
	main() 
#
