###########################################################################################################
# 
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np ##
import VPOT_conf, VPOT_1_1_VCF, VPOT_1_2_TXT #
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
info_msg1_1="VPOT: Invalid number of inputs, must have at least two :" 
info_msg1_2="VPOT: 1) output destination directory + prefix" #
info_msg1_3="VPOT: 2) input file containing VCF or TXT sample files" 
info_msg1_4="VPOT: 3) prioritisation parameters file" 
info_msg2="VPOT: Parameter file supplied" #
info_msg3_1="VPOT: A parameter file was created for these inputs files, please review and update :" #
info_msg3_2="VPOT: File is located at : " #
info_msg4="VPOT ERROR: At least one predictor or Exonic_annotation must be supplied in the parameter file. " #
#
input_type_VCF=True #
#txt_start=-1 #
#txt_end=-1 # 
sample_loc=-1 #
sample_ID="" #
INFO_loc=-1 #
FORMAT_loc=-1 #
#GT_loc=-1 #
header_ln="" #
blank_variant_ln="" #
#
############################################################################################################
#
###########################################################################################################
def initial_setup():
	global supplied_args 
	#
#	print "initial_setup():" #
#	print suffix #
#	print sys.argv #
	supplied_args=len(sys.argv) #
#	print supplied_args #
#
#	print ("args : ", supplied_args) #
	if (supplied_args < 4 ):  # arg [0] is the python program
		print (info_msg1_1+nl+info_msg1_2+nl+info_msg1_3+nl+info_msg1_4) #
		return 1 #
	else :
		VPOT_conf.output_dir=sys.argv[2] #
#		print (VPOT_conf.output_dir) #
		VPOT_conf.input_file=sys.argv[3] #
#		print (VPOT_conf.input_file) #
		if (supplied_args == 5 ):  # have supplied parameter file
			VPOT_conf.parameter_file=sys.argv[4] #
#			print ("user supplied") #
		else : ##		
			VPOT_conf.parameter_file=VPOT_conf.output_dir+"default_prioritisation_parameters.txt" #
		VPOT_conf.working_file1=VPOT_conf.output_dir+"working_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.working_file2=VPOT_conf.output_dir+"working_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.full_file1=VPOT_conf.output_dir+"full_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.full_file2=VPOT_conf.output_dir+"full_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.sort_file1=VPOT_conf.output_dir+"sort_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.sort_file2=VPOT_conf.output_dir+"sort_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.final_output_file=VPOT_conf.output_dir+"final_output_file_"+suffix+".txt" #
	
	return 0 #
#
###########################################################################################################
#
###########################################################################################################
def create_annotation_parameter(first_inputfn): #
#
#	print "create_annotation_parameter(input_file,parameter_file): " #
	annotation_header1=VPOT_conf.PF+tab+"Population_filter"+tab+"Value"+nl #
	annotation_header2=VPOT_conf.PD+tab+"Predictors"+tab+"Type"+tab+"Pred_val"+tab+"VPOT_Value"+tab+"Pred_val"+tab+"VPOT_Value"+tab+"Pred_val"+tab+"VPOT_Value"+nl # 
	annotation_header3=VPOT_conf.VT+tab+"Variant_annotation"+tab+"Exception_variant_types"+tab+"Value"+nl # 
	annotation_header4=VPOT_conf.GN+tab+"Gene Symbol"+nl # 
	annotation_header4_1=VPOT_conf.GN+tab+"Gene.refGene"+nl # 
	annotation_header5=VPOT_conf.QC+tab+"Quality Control"+tab+"Value"+nl # 
	annotation_header5_1=VPOT_conf.QC+tab+"Coverage"+tab+"0"+nl # 
	annotation_header5_2=VPOT_conf.QC+tab+"Hete_Balance"+tab+"0"+nl # 
	annotation_header5_3=VPOT_conf.QC+tab+"Genotype_Quality"+tab+"0"+nl # 
	annotation_header6=VPOT_conf.VS+tab+"Variant Score Threshold"+tab+"Value"+nl # 
	annotation_header6_1=VPOT_conf.VS+tab+"Score"+tab+"0"+nl # 
#	annotation_header6_2=VPOT_conf.VS+tab+"Percentage"+tab+"0"+nl # 
##
	#
#	print "1st file of input: ",first_inputfn #
	if (input_type_VCF) : # VCF input
		VPOT_1_1_VCF.parameters(first_inputfn) #
	else : # txt input
		VPOT_1_2_TXT.parameters(first_inputfn) #
#
#	print VPOT_conf.pop_array #
#	print VPOT_conf.pred_array #
	#
#
	VPOT_conf.startlen=len(VPOT_conf.pred_array[0]) #
#	print startlen #
	if (input_type_VCF) : #
		VPOT_1_1_VCF.setup_default_pred_values(first_inputfn) # now setup data from the 1st sample file
	else : # VCF input
		VPOT_1_2_TXT.setup_default_pred_values(first_inputfn) # now setup data from the 1st sample file
#		#
	param_file = open(VPOT_conf.parameter_file,'w',encoding="utf-8") #
#	
	param_file.write(annotation_header1) #
	for MAF in VPOT_conf.pop_array : 
		out1="\t".join(MAF)+nl #
		param_file.write(out1) #
	#
	param_file.write(annotation_header2) #
	for PRED in VPOT_conf.pred_array : 
		out1="\t".join(map(str,PRED))+nl #
		param_file.write(out1) #
	#
	#
	param_file.write(annotation_header3) #
#	print startlen 
	for PRED in VPOT_conf.pred_array : 
		if any( s in PRED[1] for s in VPOT_conf.Exonic ) : # when filtering for QC value 
			aa=len(PRED) 
			k=VPOT_conf.startlen # point to 1st option slot 
#			print aa,"/",k,"/",PRED 
			while (k < aa) :
				out1=VPOT_conf.VT+tab+PRED[1]+tab+PRED[k]+tab+"50"+nl #
				param_file.write(out1) #
				k+=1 # move to pred_array slot
#
	#
	param_file.write(annotation_header4) #
	param_file.write(annotation_header4_1) #
	param_file.write(annotation_header5) #
	param_file.write(annotation_header5_1) #
	param_file.write(annotation_header5_2) #
	param_file.write(annotation_header5_3) #
	param_file.write(annotation_header6) #
	param_file.write(annotation_header6_1) #
#	param_file.write(annotation_header6_2) #
#
	print (info_msg3_1) #
	print (info_msg3_2, VPOT_conf.parameter_file) #
#
	param_file.close() #
###########################################################################################################
#
###########################################################################################################
def read_parameter_file(): #
#
##
#	print "read_parameter_file():" #
	print (info_msg2) # parameter file supplied then use it
#	print VPOT_conf.parameter_file #
	#
	with open(VPOT_conf.parameter_file,'r',encoding="utf-8") as parm_fn : #
		for line in parm_fn: #
			this_line=re.split('\t|\n|\r|=|;',line) #
#			print line 
			if (this_line[0] == VPOT_conf.PF ): # is this a population filter line 
				if (this_line[1] != "Population_filter" ): # if not the header line - then 
						VPOT_conf.PF_array.append(this_line) 
#
			if (this_line[0] == VPOT_conf.PD ): # is this a predictor line 
				if (this_line[1] != "Predictors" ): # if not the header line - then 
						VPOT_conf.PD_array.append(this_line) 
#
			if (this_line[0] == VPOT_conf.VT ): # is this an Exonic_annotation line 
				if (this_line[1] != "Variant_annotation" ): # if not the header line - then 
						VPOT_conf.VT_array.append(this_line) 
#
			if (this_line[0] == VPOT_conf.GN ): # is this a Gene reference line 
				if (this_line[1] != "Gene Symbol" ): # if not the header line - then 
						VPOT_conf.GN_value=this_line[1] # save the gene reference field 
	#
			if (this_line[0] == VPOT_conf.QC ): # is this a Quality Control reference line 
				if (this_line[1] != "Quality Control" ): # if not the header line - then 
					if (this_line[1] == "Coverage" ): # if coverage value - then 
						VPOT_conf.Maxcoverage=int(this_line[2]) # save it 
					if (this_line[1] == "Hete_Balance" ): # if Hete_Balance - then 
						VPOT_conf.Hete_Balance=int(this_line[2]) # save it 
					if (this_line[1] == "Genotype_Quality" ): # if Genotype_Quality - then 
						VPOT_conf.Genotype_Quality=int(this_line[2]) # save it 
	#
			if (this_line[0] == VPOT_conf.VS ): # is this a Variant Score Threshold reference line 
				if (this_line[1] != "Variant Score Threshold" ): # if not the header line - then 
					if (this_line[1] == "Score" ): # if threshold value - then 
						VPOT_conf.VariantScoreThreshold=int(this_line[2]) # save it 
					if (this_line[1] == "Percentage" ): # if threshold value - then 
						VPOT_conf.VariantPercentageThreshold=int(this_line[2]) # save it 
	#

###########################################################################################################
#
###########################################################################################################
def create_final_output_file(): #
#
	print ("final_output_file: ",VPOT_conf.final_output_file) #
#
	with open(VPOT_conf.final_output_file,'w',encoding="utf-8") as ffile : # 
		outline1 = "Ranking"+tab+"Priority_Score"+tab+VPOT_conf.header_ln+nl #
		ffile.write(outline1) # write new header line with new sample ID to the new combined samples output file
#
	tmp_final_output_file=VPOT_conf.output_dir+"temp_final_output_file_"+suffix+"_tmp.txt" #
#
	COMMAND="sort -k1,1nr -k2,3V -k4,6 "+VPOT_conf.working_file1+" > "+tmp_final_output_file #  
	subprocess.call(COMMAND, shell=True) #
#
	priority_score=0 # initialise score 
	with open(tmp_final_output_file,'r',encoding="utf-8") as variants_file, open(VPOT_conf.final_output_file,'a',encoding="utf-8") as score_file : # 
		for line1 in variants_file: # work each line of new sample vcf file 
			line_parts=re.split('\t|\n|\r',line1) # split the variant up
#			print "line part 0 : ",line_parts #
#			print "priority score ",priority_score #
			if ( float(priority_score) == 0) : # possible 1st variants 
#				print "here 1",priority_score #
				if ( float(line_parts[0]) == 0 ) : # variant line also 0
#					print "here 2",line_parts[0] #
					score=0 #
				else : # set the priority 
					priority_score=line_parts[0] #
					score=float(line_parts[0])/float(priority_score) #
			else : # have a priority to work with
				score=float(line_parts[0])/float(priority_score) #
			chgln=str(round(score,2))+tab+line1 #
#
			if ( int(score*100) >= int(VPOT_conf.VariantPercentageThreshold) ) : # check priority score,is it larger than threshold
				score_file.write(chgln) #
#
#
#	COMMAND="sed -i 's/\\\\x3b/,/g' "+VPOT_conf.final_output_file #  
#	subprocess.call(COMMAND, shell=True) #
	#
#
	COMMAND="sed 's/\\\\x3b/,/g' "+VPOT_conf.final_output_file+" > "+VPOT_conf.full_file1 #  
	subprocess.call(COMMAND, shell=True) #
	#
	copyfile(VPOT_conf.full_file1,VPOT_conf.final_output_file) # copy the new full file

#
			
##
###########################################################################################################
#
###########################################################################################################
def main(): #
##
#	global header_fn 
	global input_type_VCF 
#
	print ("VPOT Prioritisation - Main") #
	#
#
	#
	VPOT_conf.init() #
#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
	#
	#
	print ("Samples input files : ",VPOT_conf.input_file) #
	with open(VPOT_conf.input_file,'r',encoding="utf-8") as input: # 
		file1=input.readline() #
		this_file=re.split('\t|\n|\r',file1) # split into file location and sample id
		filename=this_file[0].split('.') # split the file name to determine if VCF or TXT file 
		fn_len=len(filename) #
		#
		print (filename[fn_len-1]) #
		if ((filename[fn_len-1] == "txt") or (filename[fn_len-1] == "tsv")) : #
			input_type_VCF=False #
	#
	input.close() #
#
	if (supplied_args != 5): # no parameter file, so build it 
		create_annotation_parameter(this_file[0]) #
	else :
#		print info_msg2 # parameter file supplied then use it 
		read_parameter_file() #
		print ("QC MaxCOverage : ",VPOT_conf.Maxcoverage) #
		print ("QC Hete_balance : ",VPOT_conf.Hete_Balance) #
		print ("QC Genotype_Quality : ",VPOT_conf.Genotype_Quality) #
		print ("VS Score Threshold : ",VPOT_conf.VariantScoreThreshold) #
#		print ("VS Percentage Threshold : ",VPOT_conf.VariantPercentageThreshold) #
		
		if (input_type_VCF) : # VCF input
			error=VPOT_1_1_VCF.read_variant_source_file() #
		else : # txt input
			error=VPOT_1_2_TXT.read_variant_source_file() #
#
# the variant file from all the input VCF files have been created. 
# Now calculate the priority score for each of these variants 
#		print "error : ",error #
		if ( error != 1 ): #  no error 
			if ((len(VPOT_conf.PD_array) > 0 ) or (len(VPOT_conf.VT_array) > 0)) : # there is predictors/variant scoring  
				if (input_type_VCF) : # VCF input
					VPOT_1_1_VCF.score_the_variants() #
				else : # txt input
					VPOT_1_2_TXT.score_the_variants() #
				create_final_output_file() #V2
			else : #
				print (info_msg4) #
	#
##V2			create_final_output_file() #
#
###########################################################################################################
