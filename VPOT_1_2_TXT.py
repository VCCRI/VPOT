###########################################################################################################
# 
###########################################################################################################
#
import logging, sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf #
from shutil import copyfile #
#from __main__ import *
tab='\t' # 
nl='\n' #
#
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG) # to show debug message remove '#' at the beginning of line
#
###########################################################################################################
#
###########################################################################################################
def parameters(input_file):
#
#	global pop_array 
#	global pred_array 
#	global txt_start 
#	global txt_end 
#
	print ("VPOT_txt.parameters: ") #
	#
	with open(input_file,'r',encoding="utf-8") as first_fn : #
		predictors=False #
		line=first_fn.readline() # for txt file the first line contains all the annotation elements header
		this_line=re.split('\t',line) #
		len_line=len(this_line) #
		for j in range(len(this_line)): #
#			print this_line[j] #
			if (this_line[j] =="ALT") : # Alt column is #4
				predictors=True #
				VPOT_conf.txt_start=j+1 #
			elif (this_line[j] == "FORMAT") : # ok to stop checking for predictors 
#			elif ((this_line[j] == "FORMAT") and predictors ): # ok to stop checking for predictors 
				predictors=False #
				VPOT_conf.txt_end=j #
#				print "predictors = FALSE"
#
			if (predictors and this_line[j] !="Alt") : # true
				if any( s in this_line[j] for s in VPOT_conf.Population ) : # when filtering for QC value 
					VPOT_conf.pop_array.append([VPOT_conf.PF,this_line[j],VPOT_conf.MAFval]) 
				VPOT_conf.pred_array.append([VPOT_conf.PD,this_line[j],"N","","0","","1","","2"]) 
##
###########################################################################################################
#
###########################################################################################################
def setup_default_pred_values(file1): #
#
#	global pred_array 
#	global txt_start 
#	global txt_end 
#
	print ("VPOT_txt.setup_default_pred_values: ") #
##
	with open(file1,'r',encoding="utf-8") as first_fn : #
		for line in first_fn: #
			this_line=re.split('\t',line) #
#			print line 
			len_line=len(this_line) #
#			print this_line #
#			if (this_line[0] !="Chr") : # skip header line
			if (this_line[0] !="#CHROM") : # skip header line
				try: #
					for j in range(VPOT_conf.txt_start,VPOT_conf.txt_end): # keep values - start with after Alt column and ends with FORMAT 
						pred_array_pos=j-VPOT_conf.txt_start     # current VPOT_conf.pred_array position for value 
#					print "VPOT_conf.pred_array_pos :",VPOT_conf.pred_array_pos     # current VPOT_conf.pred_array position for value 
#					print "this line :",this_line[j]     # current VPOT_conf.pred_array position for value 
						if ((this_line[j] != ".") and (this_line[j] != "-") and (this_line[j] != "-999")): 
#						if (VPOT_conf.is_number(this_line[j])): # numeric
							if (VPOT_conf.is_number(this_line[j]) and (VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_type] != "A" )): # numeric
								this_numeric=True
								if ((VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_low] == "") or (float(VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_low]) > float(this_line[j]))): 
									VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_low] = this_line[j] 
								if ((VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_high] == "") or (float(VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_high]) < float(this_line[j]))): 
									VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_high] = this_line[j] 
							else : # if not a value that can be expressed as a floating point then is alpha
#							print "except", this_line[i] 
								this_numeric=False
								if ((VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_low] == "") or (VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_low] > this_line[j])): 
									VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_low] = this_line[j] 
								if ((VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_high] == "") or (VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_high] < this_line[j])): 
									VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_high] = this_line[j] 
								VPOT_conf.pred_array[pred_array_pos][VPOT_conf.PD_type] = "A" 
								aa=len(VPOT_conf.pred_array[pred_array_pos]) 
#							print aa #
								if (aa < VPOT_conf.Maxval): #still can add some more alpha options
									k=VPOT_conf.startlen # point to 1st option slot 
#								print k, aa #
									while (k <= aa) :
										if (k == aa): # completed search and not found	
											VPOT_conf.pred_array[pred_array_pos].append(this_line[j]) # then add it
#										print "add- ",this_line[i+1],"/",this_line[i] 
										else : # keep searching
											if (VPOT_conf.pred_array[pred_array_pos][k] == this_line[j] ) : # current value already added 
												break #
										k+=1 # move to pred_array slot
#								print pred_array[pred_array_pos] #
				except: # debug messages
					print("Error occurred at line :", line) #
					print("pred value in line :", this_line[j],":", this_numeric) #
#				print("pred_array when error occurred :", VPOT_conf.pred_array) #
#				print("pred_array index when error occurred :", j) #
					print("pred_array when error occurred :", VPOT_conf.pred_array[pred_array_pos]) #
					sys.exit(1) #
#
#	print pred_array #
	#
# determine a middle value for the numeric predictors 
	for j in range(len(VPOT_conf.pred_array)): #
#		print pred_array[j] 
		if (VPOT_conf.pred_array[j][1] == "N") : # this predictor is numeric 
			if (VPOT_conf.is_number(VPOT_conf.pred_array[j][VPOT_conf.PD_high]) and VPOT_conf.is_number(VPOT_conf.pred_array[j][VPOT_conf.PD_low])): #
				mid_value=float((float(VPOT_conf.pred_array[j][VPOT_conf.PD_high])+float(VPOT_conf.pred_array[j][VPOT_conf.PD_low]))/2) #
#			print "mid-",mid_value #
				VPOT_conf.pred_array[j][VPOT_conf.PD_mid] = mid_value 
#
#	print VPOT_conf.pred_array 

###########################################################################################################
#
###########################################################################################################
def read_variant_source_file(): #
#
##
#	print "read_variant_source_file_for_txt(): " #
#	print info_msg2 # parameter file supplied then use it
	print (VPOT_conf.parameter_file) #
	#
	with open(VPOT_conf.input_file,'r',encoding="utf-8") as input: # 
		for line in input: # work each input vcf file 
			this_line=re.split('\t|\n|\r',line) # split into file location and sample id
#			print (this_line) #
			if ( setup_for_this_src_file(this_line) != 0 ): # check the input file 
				print ("Issue with input file 1-: ", line) # 
				return 1 # error return
			else : #
				print ("processing input file- : ", this_line) # 
#				print ("checking values: ",VPOT_conf.sample_loc,VPOT_conf.FORMAT_loc,VPOT_conf.sample_coverage_loc[0]) #
				if ((VPOT_conf.sample_loc >=0) and (VPOT_conf.FORMAT_loc >=0) and (VPOT_conf.sample_coverage_loc[0] >=0)): #
#					print "OK" #
					VPOT_conf.header_ln = VPOT_conf.header_ln+tab+VPOT_conf.sample_ID # add sample ID to header line
#					print "header- ",VPOT_conf.header_ln #
#					outline1 = VPOT_conf.header_ln+nl #
					VPOT_conf.blank_variant_ln = VPOT_conf.blank_variant_ln+tab+"0" # this is a default variant line for samples that do not have this variant
#					print "blank_variant_ln- ",VPOT_conf.blank_variant_ln #
#					print sample_coverage, "/", sample_coverage_loc #
					work_this_src_file(this_line)  #
					if (os.path.isfile(VPOT_conf.full_file2) is False ) or (os.stat(VPOT_conf.full_file2).st_size == 0): # are there any variants in the final file? 
						VPOT_conf.update_existing_variants(VPOT_conf.working_file1,"1") # No - then create it based on sample file 
					elif (os.stat(VPOT_conf.working_file1).st_size == 0): # there any no variants for this sample 
						copyfile(VPOT_conf.full_file2,VPOT_conf.working_file1) # copy the new full file
						VPOT_conf.update_existing_variants(VPOT_conf.working_file1,"0") # then for existing variants, just add 0 to the end for this sample 
					else : # yes, then incorporate them 
						VPOT_conf.incorporate_this_src_into_full_file() #
					copyfile(VPOT_conf.full_file2,VPOT_conf.full_file1) # copy the new full file
					
					
				else : #
					print ("Issue with input file 2-: ", line) #
					return 1 # error return 
###########################################################################################################
#
###########################################################################################################
def setup_for_this_src_file(file_line): #
#
#	global GT_loc #
	
#
#	print ("setup_for_this_src_file for txt: ", file_line )#
#	print info_msg2 # parameter file supplied then use it
	VPOT_conf.sample_loc=-1 #
	VPOT_conf.sample_ID=""
#	VPOT_conf.INFO_loc=-1 #
	VPOT_conf.FORMAT_loc=-1 #
	VPOT_conf.sample_coverage_loc = [-1,-1,-1,-1,-1,-1] # location of VCF format codes for sample 
	#
	with open(file_line[0],'r',encoding="utf-8") as source_vcf : # open the sample input file
		for src_line in source_vcf: # work each line of source file 
			src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
#			print "src_line : ",src_line1 #
			if ("#CHROM" in src_line1[0]): # find sample location - using the header line
				for i, content in enumerate(src_line1): # return the value and index number of each item in the line array 
#								print "content-",content,"/",i				#
					if (content == file_line[1]) : # when filtering for sample ID 
						VPOT_conf.sample_loc=i #save sample location
						VPOT_conf.sample_ID=file_line[1] # save sample ID
#						print "SAMPLE: ",VPOT_conf.sample_loc,"/",VPOT_conf.sample_ID #
					elif (content == "FORMAT") : # 
						VPOT_conf.FORMAT_loc=i #save FORMAT location
#						print "FORMAT_loc: ",FORMAT_loc #
				#	
				if ("#CHROM" not in VPOT_conf.header_ln): # is this first time
#					print ("setup header") #
					VPOT_conf.header_ln='\t'.join(src_line1[:VPOT_conf.FORMAT_loc+1])+tab+VPOT_conf.GENE_NAME # save heading up to INFO (using INFO_loc +1 to denote the stop field)
#					VPOT_conf.header_ln=VPOT_conf.header_ln+tab+VPOT_conf.GENE_NAME # 
#					VPOT_conf.blank_variant_ln='\t'.join(src_line1[:VPOT_conf.FORMAT_loc+1]) # setup dummy variant line
					VPOT_conf.blank_variant_ln=VPOT_conf.header_ln # save heading up to INFO (using INFO_loc +1 to denote the stop field)
#					print "header_ln : ",VPOT_conf.header_ln #
			else : # variants lines 
#				print ("var line : ",src_line1) #
				FORMAT1=re.split(':',src_line1[VPOT_conf.FORMAT_loc]) # split the variant line's FORMAT field
#				print ("FORMAT : ",src_line1[VPOT_conf.FORMAT_loc]) # split the variant line's FORMAT field
				for i, content in enumerate(FORMAT1): # return the value and index number of each item in the line array 
#					print ("content-",content,"/",i)				#
					for j in range(len(VPOT_conf.sample_coverage)): #
#						print (VPOT_conf.sample_coverage[j],"/",content) #
						if (content == VPOT_conf.sample_coverage[j]) : 	# look for the coverage FORMAT fields 
#							print (i,"/",j,"/",VPOT_conf.sample_coverage_loc)
							VPOT_conf.sample_coverage_loc[j]=i 			# save its location	
							break #
#				print (i,"/",j)
#				print ("coverage : ",VPOT_conf.sample_coverage_loc) #
				source_vcf.close() # finish with source file 
				return 0 # have setup all location values - ok to go back
#						print INFO1 #
#						print FORMAT1 #
#						print SAMPLE1 #

#ed	source_vcf.close() # finish with source vcf file 
	return 1 # not a clean exit
#
###########################################################################################################
#
###########################################################################################################
def work_this_src_file(file_line): #
#
##
#	print "work_this_src_file(file_line): " #
#	print working_file1 #
	wrkf1=open(VPOT_conf.working_file1,'w',encoding="utf-8") # 
	with open(file_line[0],'r',encoding="utf-8") as source_vcf : #
		for src_line in source_vcf: # work each line of source vcf file 
			src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
			if ("#CHROM" not in src_line1[0]): # skip the header lines
#				print (src_line1) #
# variants lines 
				SAMPLE1=re.split(':',src_line1[VPOT_conf.sample_loc]) # split the sample's FORMAT fields 
#				print ("MaxCOverage : ",VPOT_conf.Maxcoverage) #
#				print ("Hete_balance : ",VPOT_conf.Hete_Balance) #
#				print "coverage_loc : ",VPOT_conf.sample_coverage_loc #
#				print SAMPLE1 #
#				print VPOT_conf.sample_coverage_loc #
# have a coverage check on the sample 
#				if (((VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val] == -1) or ((SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] != ".") and (int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]]) >= VPOT_conf.Maxcoverage))) and #  NR
#					((VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val] == -1) or ((SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]] != ".") and (int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]]) >= VPOT_conf.Maxcoverage)))) : # DP
#				if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "./.") : # a valid genotype 
				if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "./.") and (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "0/.") and (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "./0") and (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "0/0")  : # a valid genotype 
					if (VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val] != -1 ) : #this sample have a coverage depth
						Sample_coverage=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]]) # save DP value
						Alt_reads=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]])/2 # save DP value
#						print ("DP") #
					if (VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val] != -1 ) and (VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val] != -1 ) : #this sample have a coverage depth from NR and NV
#						Sample_coverage=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]])+int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]]) # save DP value
						Sample_coverage=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]]) # save DP value
						Alt_reads=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]]) # save DP value
#						print ("NR+NV") #
#					print ("TOT: ",str(Sample_coverage)) # 
#					print ("ALT_READ: ",str(Alt_reads)) # 
#					print ("pass") #
#					if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]] == ".") : # no DP_val 
#						SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]] = "0"  # set it as zero 
#					if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] == ".") : # no DP_val 
#						SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] = "0"  # set it as zero 
#					if (((VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val] == -1) or 
#						((VPOT_conf.is_number(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]])) and (int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]]) >= int(VPOT_conf.Maxcoverage)))) and  #  NR
#						((VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val] == -1) or 
#						((VPOT_conf.is_number(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]])) and (int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]]) >= int(VPOT_conf.Maxcoverage))))) : # DP
#
					VPOT_conf.QC_PASS=False #
					if ((Sample_coverage >= int(VPOT_conf.Maxcoverage)) and (int((Alt_reads/Sample_coverage)*100) >= int(VPOT_conf.Hete_Balance))) : # Pas QC for coverage and balance
						VPOT_conf.QC_PASS=True # Yes
#					print ("QC_PASS",VPOT_conf.QC_PASS) #
#						print ("add") #
					GT_values=re.split('/',SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]]) # get the genotype fields
# 	 				print GT_values #
					for j in range(len(GT_values)) : #
#					print GT_values[j] #
						if ( GT_values[j] not in VPOT_conf.Non_alt_GT_types ) : # when filtering for QC value 
#							print ("keep this variant1") #
							check_this_variant(src_line, wrkf1) #
							break # get out of for loop (GT_values)

	wrkf1.close() # finish with the output file 
#
#	print "sort unique" #
##	COMMAND="sort -V -u -k1,5 "+VPOT_conf.working_file1+" > "+VPOT_conf.sort_file1 #  
	COMMAND="sort -V -u "+VPOT_conf.working_file1+" > "+VPOT_conf.sort_file1 #  
	subprocess.call(COMMAND, shell=True) #
	copyfile(VPOT_conf.sort_file1,VPOT_conf.working_file1) # copy back

#	source_vcf.close() # finish with source vcf file 
#
###########################################################################################################
#
###########################################################################################################
def check_this_variant(src_line, wrkf1):  #
#
#	print "check_this_variant(src_line, wrkf1):  ",src_line #
#
	src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
	SAMPLE1=re.split(':',src_line1[VPOT_conf.sample_loc]) # split the sample's FORMAT fields 
	GENOTYPE1=re.split('/',SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]]) # split the sample's GENOTYPE fields 
	header1=re.split('\t|\n|\r',VPOT_conf.header_ln) # split into file location and sample id
#	print "header : ",VPOT_conf.header_ln #
#	print (src_line1) #
	if (population_frequency(src_line1,header1) == 0 ): # check if variant with in filter
#		print ("saving this variant1") #
		VPOT_conf.gene_ref="NONE" #
		find_gene_ref(src_line1,header1)    # check if variant found != means yes 
#		if ( gene != 0 ):    # check if variant found != means yes 
#			VPOT_conf.gene_ref=src_line1[gene] # now add in the gene ref 
#
		if (not VPOT_conf.QC_PASS) : # sample did not pass QC  
			gtype="." #
#		elif (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] == "1/1") : # a homozygous alt genotype 
		elif (GENOTYPE1[0] == GENOTYPE1[1]) : # a homozygous alt genotype 
			gtype="2" #
		else : #
			gtype="1" 
#
		outline='\t'.join(src_line1[:VPOT_conf.FORMAT_loc+1])+tab+VPOT_conf.gene_ref+tab+gtype+nl #
		wrkf1.write(outline) #
#	else : #
#		print "did not pass PF" #
#
###########################################################################################################
#
###########################################################################################################
def population_frequency(info_ln,header1): #
#
#
#	print ("population_frequency(info_ln): ") #
#	print "header : ",header1 #
#	print info_ln[1] #
	#
	val=0 #
#	INFO1=re.split('\t',info_ln) # split into file location and sample id
#	print (VPOT_conf.PF_array) #
	for j in range(len(VPOT_conf.PF_array)): #
#		print (VPOT_conf.PF_array[j][1]) # move to pred_array slot
		for i, content in enumerate(header1): # return the value and index number of each item in the line array 
#			print ("content-",content,"/",i,"/",VPOT_conf.PF_array[j][1]) 				#
			if (content == VPOT_conf.PF_array[j][1]) : # the population freq annotation we want? 
				if ( info_ln[i] != ".") :	 # if numeric check, else ok as most likely "." to state no annotation  
					temp_val=info_ln[i]
				else : # not number, so set a default as no annotation in population_frequency check would mean novel	
					temp_val=-999
				if ( VPOT_conf.is_number(temp_val) ) : # numeric
				#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
#					print ("temp_val yes :",temp_val) #
#					print (VPOT_conf.PF_array[j][3]) #
					if ( VPOT_conf.PF_array[j][3] == "LE" ) :	 # lower or each to PF limit 
#						print (info_ln[i],"/",VPOT_conf.PF_array[j][2]) #
##				if ( float(info_ln[i]) > float(VPOT_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
						if ( float(temp_val) > float(VPOT_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
#						print ("do not want : ",info_ln[i],"/",VPOT_conf.PF_array[j][2]) #
							val=1 # do not want this variant 
							break #
					else : # GT  look for greater or equal to PF limit	
						if ( float(temp_val) < float(VPOT_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
							val=1 # do not want this variant 
							break #
#
	return val #
#	
###########################################################################################################
#
###########################################################################################################
def find_gene_ref(info_ln,header1): #
#
#
#	print "find_gene_ref(info_ln,header1): " #
#	print "header : ",header1 #
#	print info_ln[1] #
	#
	val=0 #
	for i, content in enumerate(header1): # return the value and index number of each item in the line array 
#		print "content-",content,"/",i,"/",VPOT_conf.PF_array[j][1] 				#
		if (content == VPOT_conf.GN_value ) : # the gene ref annotation we want? 
			VPOT_conf.gene_ref=info_ln[i] # save gene name 
			break #
#	
###########################################################################################################
#
###########################################################################################################
def score_the_variants(): #
##
#	print "score_the_variants(): " #
	#
	header1=re.split('\t|\n|\r',VPOT_conf.header_ln) # split into file location and sample id
	with open(VPOT_conf.full_file1,'r',encoding="utf-8") as variants_file, open(VPOT_conf.working_file1,'w',encoding="utf-8") as score_file : # 
		for line1 in variants_file: # work each line of new sample vcf file 
			priority_score=0 # initialise score 
			line_parts=re.split('\t|\n|\r',line1) # split the variant up
#			print ("line parts : ",line_parts) #
#			print ("line part 0 : ",line_parts[0]) #
			if ("#CHROM" != line_parts[0]): #
#				print src_line1 #
				logging.debug('variant line : %s', str(line1))
#				
				if (len(VPOT_conf.PD_array) > 0 ) : #  
					priority_score=prioritise_variants_by_predictors(line_parts,header1) # get priority score
				#
				if (len(VPOT_conf.VT_array) > 0 ) : #   
#					print ("check VT") #
					priority_score=priority_score+prioritise_variants_by_VT_types(line_parts,header1) # 
				#
				if ( float(priority_score) <= 0 ) : # check priority score
					priority_score=0 #
#
				if ( int(priority_score) >= int(VPOT_conf.VariantScoreThreshold) ) : # check priority score,is it larger than threshold
					outline=str(priority_score)+tab+line1 # yes - then output it
					score_file.write(outline) #
#				print outline #
#				score_file.write(outline) #
			else : # save the header line	
				outline = "Priority_score"+tab+line1 #
				score_file.write(outline) #
#
###########################################################################################################
#
###########################################################################################################
def prioritise_variants_by_predictors(INFO_details,header1): #
#
#	global PD_array 
#
#	print "prioritise_variants_by_predictors(INFO_details): " #
#
	type=2 # array location that contains the value for the type of predictor (alpha or numeric)
	a_value=a_score=0 # array location variable
#
	val=wrkval=pdval=0 # this is the predictor priority score
#	print INFO1 #
	for j in range(len(VPOT_conf.PD_array)): #
#		print (VPOT_conf.PD_array[j][1]) # move to pred_array slot
		pdval=0 #
		for i, content in enumerate(header1): # return the value and index number of each item in the line array 
#			print "content-",content,"/",i				#
#			print content #
			if (content == VPOT_conf.PD_array[j][1]) : # the predictor annotation we want? 
				limit=len(VPOT_conf.PD_array[j]) # 
#				print "PD_array: ", limit, PD_array[j] #
#				print a_score, "/", a_value #
				a_value=type+1 #
				a_score=a_value+1 #
#				print ("INFO_details : ",INFO_details[i]) #
				pdval=0 #
				predval=re.split('\||&',INFO_details[i]) # split into individual values using delimiter | and & 
				logging.debug('PD info : %s - %s', content, predval)
#				print ("predval : ",predval) #
				for k in range(len(predval)): #
					b_value=a_value #
					b_score=a_score #
					wrkval=0 #
					if (VPOT_conf.PD_array[j][type] == "A") : # character field 
						while (b_value < len(VPOT_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
							if ( predval[k] == VPOT_conf.PD_array[j][b_value] ) :		# is this the variant type we are looking for   
								wrkval=int(VPOT_conf.PD_array[j][b_score])	 		# yes - return the score for the variant 
								break 										# finish 
							b_value+=2	 # bump to next value option   
							b_score+=2	 #   
						if (wrkval > pdval ) : #check for higher predictor score
							pdval=wrkval	
#						while (a_value < len(VPOT_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#							if ( INFO_details[i] == VPOT_conf.PD_array[j][a_value] ) :		# is this the variant type we are looking for   
#								val=val+int(VPOT_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
#								break 										# finish 
#							a_value+=2	 # bump to next value option   
#							a_score+=2	 #   
					else : # numeric field 
#						if ( predval[k] != ".") :	 # if annotation then continue   
						if ( VPOT_conf.is_number(predval[k]) ) : # numeric
							while (b_value < len(VPOT_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
								if ( float(predval[k]) < float(VPOT_conf.PD_array[j][b_value])) :	 # is this the variant type we are looking for   
									wrkval=int(VPOT_conf.PD_array[j][b_score])	 		# yes - return the score for the variant 
									break 										# finish 
								if ( b_value+2 >= len(VPOT_conf.PD_array[j])-1) : # end of this PD limits  
#								print ("checking last value") #
									if ( float(predval[k]) > float(VPOT_conf.PD_array[j][b_value])) :	 # check the last value, which is a greater then check   
										wrkval=int(VPOT_conf.PD_array[j][b_score])	 		# yes - return the score for the variant 
								b_value+=2	 # bump to next value option   
								b_score+=2	 #   
							if (wrkval > pdval ) : #check for higher predictor score
								pdval=wrkval	
#						if ( INFO_details[i] != ".") :	 # if annotation then continue   
#							while (a_value < len(VPOT_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#								if ( float(INFO_details[i]) < float(VPOT_conf.PD_array[j][a_value])) :	 # is this the variant type we are looking for   
#									val=val+int(VPOT_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
#									break 										# finish 
#								if ( a_value+2 >= len(VPOT_conf.PD_array[j])-1) : # end of this PD limits  
#								print ("checking last value") #
#									if ( float(INFO_details[i]) > float(VPOT_conf.PD_array[j][a_value])) :	 # check the last value, which is a greater then check   
#										val=val+int(VPOT_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
#								a_value+=2	 # bump to next value option   
#								a_score+=2	 #   
#				print ("variant's PD score : ",pdval) #
				val=val+pdval	 		# add to overall variant score 
				break # done for this predictor annotation 
#
	return val #
#	
##
#	print info_msg2 # parameter file supplied then use it
###########################################################################################################
#
###########################################################################################################
def prioritise_variants_by_VT_types(INFO_details,header1): #
#
#	global VT_array 
	variant_INFO_array=[]
	variant_INFO_val_array=[]
#
	val=wrkval=pdval=-999 # this is the predictor priority score
	VT_already_stored=False# 
	VT_found_for_variant=False #
	VT_k_found=False #
	VT_k_has_annotation=False #
	VT_has_annotation_value_in_variant=False #
#	print ("Variant INFO :",INFO_details) #
	VT_array_size=len(VPOT_conf.VT_array)
#	print ("number of VT lines in PPF :",VT_array_size) #
#
#   loop thru the PPF VT lines
	for j in range(len(VPOT_conf.VT_array)): # move thru the VT lines in PPF
#		print ("VT array :", VPOT_conf.VT_array[j]) # move to pred_array slot
#
# 		loop thru the variant's INFO details, working one VT PPF line at a time
		for i, content in enumerate(header1): # return the value and index number of each item in the line array 
#			print ("content-",content,"/",i)				#
			if (content == VPOT_conf.VT_array[j][1]) : # the variant functional annotation we want? then save the annotation details for later use
				VT_already_stored=False# 
				if (len(variant_INFO_array) > 0) : # have I saved something already- yes
					for v in range(len(variant_INFO_array)): # move thru the saved VT INFO fields
						if (content == variant_INFO_array[v]) : # my current VT INFO is already saved? - yes, no need to save it again
							VT_already_stored=True# 
							break #
				if (not VT_already_stored) : # if not already saved - then save the INFO field name and values
						variant_INFO_array.append(content)
						variant_INFO_val_array.append(INFO_details[i])
	if (len(variant_INFO_array) < 1) : # NO VT INFO fields found
#		print ("no VT INFO field to match the VT PPF lines") #
		val=0	# then set a 0 score,
		return val # no further processing
######
#	print ("VT WORK array :",variant_INFO_array, variant_INFO_val_array,"\n") #
	logging.debug('VT WORK array : %s - %s', variant_INFO_array, variant_INFO_val_array)
#	global VT_array 
#
#####
# now work these INFO fields
	for v in range(len(variant_INFO_array)): # move thru the saved VT INFO fields
#		print ("VT array :", VPOT_conf.VT_array[j]) # move to pred_array slot
#		print ("content-",variant_INFO_array[v])				#
#		print ("INFO_details : ",variant_INFO_val_array[v]) #
		predval=re.split('\||&',variant_INFO_val_array[v]) # split into individual values using delimiter | and & 
#
		for k in range(len(predval)): # for the variant info fields
			VT_k_has_annotation=False # initalise for this set of VT annotaion
			VT_k_found=False # initalise for this set of VT annotaion
#
			for j in range(len(VPOT_conf.VT_array)): # move thru the VT lines in PPF
				if (variant_INFO_array[v] == VPOT_conf.VT_array[j][1]) : # the variant functional annotation we want? 
#					print ("compare VT : vcf annotation - ",predval[k],"/ PPF VT- ",VPOT_conf.VT_array[j])  #
					if (predval[k] != ".") : # an annotated value
						VT_k_has_annotation=True #
						if ( predval[k] == VPOT_conf.VT_array[j][2]) :	 # is this the specific variant type we are looking for   
							VT_k_found=True # found a match
#							print ("vtfound") #
							wrkval=int(VPOT_conf.VT_array[j][3])	 # yes - return the score for the variant
							if (wrkval > val ) : #check for higher predictor score
#								print ("new variant's VT score : ",wrkval) #
								val=wrkval	
							break # stop loop as this one has found a match
					else : # there is no annotation for this INFO field
						break
				# end of for loop
				#
#			print ("setting the inds - found -", VT_k_found, "anno -", VT_k_has_annotation,"\n") # This annotation was found in PPF and scored
			if (VT_k_found and VT_k_has_annotation) : # This annotation was found in PPF and scored
				VT_found_for_variant=True
			elif (VT_k_has_annotation) : # only had a valid annotation but not in PPF
				VT_has_annotation_value_in_variant=True
#
# have finished checking all the VT lines in the PPF against variant
# now determine what is the prority score to pass back for the variant.
#
#	print ("variant's before correction VT score : ",val) # this is the score based on matches to the PPF VT lines
#
#	print ("final inds - VT_found -", VT_found_for_variant, "/ ANNO exist in variant -", VT_has_annotation_value_in_variant,"\n") # This annotation was found in PPF and scored
	logging.debug('final inds - VT_found - %s / ANNO exist in variant - %s ', VT_found_for_variant, VT_has_annotation_value_in_variant)
	if (not VT_found_for_variant and VT_has_annotation_value_in_variant and (val == -999 )) : # variant's annotations were not in any of the VT PPF 
#		print ("no VT annotations matched in the PPF VT lines for the variant") #
		val=0	# then set a 0 score,
	elif (not VT_found_for_variant and not VT_has_annotation_value_in_variant ) : # no annotation found in the variant 
#		print ("no VT annotations exist in the variant") #
		val=0	# then set a 0 score,
	elif (VT_found_for_variant and (val < 0 )) : # there was a match but only for low score
#		print ("there was a match but only for low score") #
		if (VT_has_annotation_value_in_variant) : # if there were other annotated values in the variant
#			print ("but there were other VT annotation for variants not in PPF") #
			val=0	# then set a 0 score,
#		else : # no other annotated field, so return the low value	
#			print ("so return matched low value VT fields") #
#	elif (VT_found_for_variant and (val > 0 )) : # there was a match VT value from the PPF to the variant, and with a positive score 
#		print ("matched VT fields") #
#	print ("variant's final VT score : ",val) #
###
	return val #
##
#