###########################################################################################################
# 
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf #
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
#
Sample_loc=-1 #
#
info_msg1_1="VPOT samplef : Invalid number of inputs, must have three :" 
info_msg1_2="VPOT samplef : 1) output destination directory + prefix" #
info_msg1_3="VPOT samplef : 2) Input file - output from VPOT prioritisation process" 
info_msg1_4="VPOT samplef : 3) Selection Criteria list (ped format)" #
#
############################################################################################################
###########################################################################################################
#
###########################################################################################################
def initial_setup():
#
	global supplied_args 
	#
	print ("initial_setup():") #
	print (suffix) #
#	print sys.argv #
	supplied_args=len(sys.argv) #
#	print supplied_args #
#
	if (supplied_args != 5 ):  # arg [0] is the python program
		print (info_msg1_1+nl+info_msg1_2+nl+info_msg1_3+nl+info_msg1_4) #
		return 1 #
	else :
		VPOT_conf.output_dir=sys.argv[2] #
		VPOT_conf.input_file=sys.argv[3] #
		VPOT_conf.selection_list=sys.argv[4] #
		VPOT_conf.final_output_file=VPOT_conf.output_dir+"variant_filtered_output_file_"+suffix+".txt" #
		print ("output : ",VPOT_conf.final_output_file) #
	
	return 0 #
#
#
###########################################################################################################
def setup_samples(): #
#
#	print "setup_samples(): " #
#
	with open(VPOT_conf.selection_list,'r',encoding="utf-8") as select_input: # 
		for line1 in select_input: # work each line of new sample vcf file 
			this_line=re.split('\t|\n|\r',line1,7) # split into sample id and status (1/0)
			this_line[6]=0 # 
#			print (this_line) #
			if (this_line[5]=="2"): # is this affected sample
#				print ("affected",this_line) #
				this_line[5]="1" # set on
			else :  #
#				print ("unaffected",this_line) #
				this_line[5]="0" # set on
			VPOT_conf.Sample_ids.append(this_line) 
		#
#	print (VPOT_conf.Sample_ids) #
#
###########################################################################################################
#
#
###########################################################################################################
def filter_the_variants(): #
#
# input file for filtering is the output from the VPOT process. This means the variant gene name is in a column named GENE_NAME.
	global Sample_loc #
#
#	print "filter_the_variants(): " #
	#
	with open(VPOT_conf.input_file,'r',encoding="utf-8") as variants_file, open(VPOT_conf.final_output_file,'w',encoding="utf-8") as filtered_file : # 
		for line1 in variants_file: # work each line of new sample vcf file 
			write_it=False # initialise score 
			line_parts=re.split('\t|\n|\r',line1) # split the variant up
#			print "line part 0 : ",line_parts[0] #
			if ("#CHROM" != line_parts[2]): #
#				print src_line1 #
				write_it=filter_variants_by_Seg(line_parts) # check get priority score
				#
			else : # save the header line	
				write_it=True # initialise score 
				for i, content in enumerate(line_parts): # return the value and index number of the sample id item in the line array 
#					print "content-",content,"/",i				#
					for j in range(len(VPOT_conf.Sample_ids)): #
						if (content == VPOT_conf.Sample_ids[j][1]) : # look for sample id in header 
							VPOT_conf.Sample_ids[j][6]=i #save sample location
#						print "INFO_loc: ",INFO_loc #
					#	
				print (VPOT_conf.Sample_ids) #
			if (write_it): #
				filtered_file.write(line1) # write the line to final output file 
#
###########################################################################################################
#
###########################################################################################################
def filter_variants_by_Seg(INFO_details): #
#
#	print "filter_variants_by_GN(INFO_details): " #
#
	val=True #
	for j in range(len(VPOT_conf.Sample_ids)): #
		if ( VPOT_conf.Sample_ids[j][6] != 0 ) : # check if this samples is in the input file 
			if ( VPOT_conf.Sample_ids[j][5] != INFO_details[VPOT_conf.Sample_ids[j][6]] ) : # yes - check the sample's value 
				val=False # not what is needed 
				break # then get out and move to next variant
#
	return val #
#	
##
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	VPOT_conf.init() #
	#
	print ("Variant Segregation Filter - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
	#
# Now filter the input file by gene list 
	setup_samples() #
	filter_the_variants() #
#
#	print (VPOT_conf.Sample_ids) #
	#
#	clean_up() #
#
############################################################################################################