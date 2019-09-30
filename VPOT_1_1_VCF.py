###########################################################################################################
# 
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf #
import gzip #
from shutil import copyfile #
#from __main__ import *
tab='\t' # 
nl='\n' #
###########################################################################################################
#
###########################################################################################################
def parameters_p1(input_file):
#
#	global pop_array 
#	global pred_array 
#
		predictors=False #
		while True:
			line=input_file.readline() # go thru the INFO annotation section of the VCF header
#		print line 
			this_line=re.split('\t|\n|\r|=|,',line) #
			if (this_line[0] =="#CHROM") : 
				break #
			elif (this_line[2] =="ANNOVAR_DATE") : #
				predictors=True #
			elif (this_line[2] =="ALLELE_END") : #
				predictors=False #
#
			if (predictors and this_line[2] !="ANNOVAR_DATE") : # true
#				print this_line[2] #
#				print VPOT_conf.Population # when filtering for QC value 
				if any( s in this_line[2] for s in VPOT_conf.Population ) : # when filtering for QC value 
					VPOT_conf.pop_array.append([VPOT_conf.PF,this_line[2],VPOT_conf.MAFval,"LE"]) 
				else :
					VPOT_conf.pred_array.append([VPOT_conf.PD,this_line[2],"N","","0","","1","","2"]) 
#
###########################################################################################################
#
###########################################################################################################
def parameters(input_file):
#
#	global pop_array 
#	global pred_array 
#
	print ("VPOT_VCF.parameters: ") #
#
	#
	if input_file.endswith('.gz'):
		with gzip.open(input_file,'rt') as first_fn : #
			parameters_p1(first_fn) #
	else:
		with open(input_file,'r',encoding="utf-8") as first_fn : #
			parameters_p1(first_fn) #
###########################################################################################################
#
###########################################################################################################
def setup_default_pred_values_1(first_fn): #
#
#	global pred_array 
#
	for line in first_fn: #
		this_line=re.split('\t|\n|\r|=|;',line) #
#			print line 
		try: #
			for i, content in enumerate(this_line): # return the value and index number of each item in the line array 
#				print "content-",content,"/",i				#
				for j in range(len(VPOT_conf.pred_array)): #
					if (content == VPOT_conf.pred_array[j][1]) : # when filtering for QC value 
#						print VPOT_conf.pred_array[j][1], content, this_line[i+1] #
						#if ((this_line[i+1] != ".") and (this_line[i+1] != "-999")): 
						if ((this_line[i+1] != ".") and (this_line[i+1] != "-") and (this_line[i+1] != "-999")): # numeric		
							if (VPOT_conf.is_number(this_line[i+1])): # numeric
								if ((VPOT_conf.pred_array[j][VPOT_conf.PD_low] == "") or (float(VPOT_conf.pred_array[j][VPOT_conf.PD_low]) > float(this_line[i+1]))): 
									VPOT_conf.pred_array[j][VPOT_conf.PD_low] = this_line[i+1] 
								if ((VPOT_conf.pred_array[j][VPOT_conf.PD_high] == "") or (float(VPOT_conf.pred_array[j][VPOT_conf.PD_high]) < float(this_line[i+1]))): 
									VPOT_conf.pred_array[j][VPOT_conf.PD_high] = this_line[i+1] 
							else : # if not a value that can be expressed as a floating point then is alpha
#								print "except", this_line[i] 
								if ((VPOT_conf.pred_array[j][VPOT_conf.PD_low] == "") or (VPOT_conf.pred_array[j][VPOT_conf.PD_low] > this_line[i+1])): 
									VPOT_conf.pred_array[j][VPOT_conf.PD_low] = this_line[i+1] 
								if ((VPOT_conf.pred_array[j][VPOT_conf.PD_high] == "") or (VPOT_conf.pred_array[j][VPOT_conf.PD_high] < this_line[i+1])): 
									VPOT_conf.pred_array[j][VPOT_conf.PD_high] = this_line[i+1] 
								VPOT_conf.pred_array[j][VPOT_conf.PD_type] = "A" 
								aa=len(VPOT_conf.pred_array[j]) 
#								print aa #
								if (aa < VPOT_conf.Maxval): #still can add some more alpha options
									k=VPOT_conf.startlen # point to 1st option slot 
#									print k, aa #
									while (k <= aa) :
										if (k == aa): # completed search and not found	
											VPOT_conf.pred_array[j].append(this_line[i+1]) # then add it
#											print "add- ",this_line[i+1],"/",this_line[i] 
										else : # keep searching
											if (VPOT_conf.pred_array[j][k] == this_line[i+1] ) : # current value already added 
												break #
										k+=1 # move to VPOT_conf.pred_array slot
#									print VPOT_conf.pred_array[j] #
		except: # debug messages
			print("Error occurred at line :", line) #
			print("pred value in line :", this_line[i+1]) #
#			print("pred_array when error occurred :", VPOT_conf.pred_array) #
#			print("pred_array index when error occurred :", j) #
			print("pred_array when error occurred :", VPOT_conf.pred_array[j]) #
			sys.exit(1) #
#
###########################################################################################################
#
###########################################################################################################
def setup_default_pred_values(file1): #
#
#	global pred_array 
#
	print ("VPOT_VCF.setup_default_pred_values: ") #
##
	if file1.endswith('.gz'):
		with gzip.open(file1input_file,'rt') as first_fn : #
			setup_default_pred_values_1(first_fn) #
	else:
		with open(file1,'r',encoding="utf-8") as first_fn : #
			setup_default_pred_values_1(first_fn) #
#
#	print pop_array #
#	print VPOT_conf.pred_array #
	#
# determine a middle value for the numeric predictors 
	for j in range(len(VPOT_conf.pred_array)): #
#		print VPOT_conf.pred_array[j] 
		if (VPOT_conf.pred_array[j][1] == "N") : # this predictor is numeric 
			if (VPOT_conf.is_number(VPOT_conf.pred_array[j][VPOT_conf.PD_high]) and VPOT_conf.is_number(VPOT_conf.pred_array[j][VPOT_conf.PD_low])): #
				mid_value=float((float(VPOT_conf.pred_array[j][VPOT_conf.PD_high])+float(VPOT_conf.pred_array[j][VPOT_conf.PD_low]))/2) #
#			print "mid-",mid_value #
				VPOT_conf.pred_array[j][VPOT_conf.PD_mid] = mid_value 
#
#	print VPOT_conf.pred_array 

#
###########################################################################################################
#
###########################################################################################################
def read_variant_source_file(): #
#
##
#	print "read_variant_source_file_for_VCF(): " #
#	print info_msg2 # parameter file supplied then use it
	print (VPOT_conf.parameter_file) #
	#
	with open(VPOT_conf.input_file,'r',encoding="utf-8") as input: # 
		for line in input: # work each input vcf file 
			this_line=re.split('\t|\n|\r',line) # split into file location and sample id
#			print this_line #
			if ( setup_for_this_src_file(this_line) != 0 ): #
				print ("Issue with input file 1: ", line) #
				return 1 # error

			else : #
				print ("processing input file : ", this_line) # 
				if ((VPOT_conf.sample_loc >=0) and (VPOT_conf.INFO_loc >=0) and (VPOT_conf.FORMAT_loc >=0) and (VPOT_conf.sample_coverage_loc[0] >=0)): #
#					print "OK" #
					VPOT_conf.header_ln = VPOT_conf.header_ln+tab+VPOT_conf.sample_ID # add sample ID to header line
#					print "header- ",header_ln #
#					outline1 = VPOT_conf.header_ln+nl #
					VPOT_conf.blank_variant_ln = VPOT_conf.blank_variant_ln+tab+"0" #
#					print "blank_variant_ln- ",blank_variant_ln #
#					print sample_coverage, "/", sample_coverage_loc #
					work_this_src_file(this_line)  #
					if (os.path.isfile(VPOT_conf.full_file2) is False ) or (os.stat(VPOT_conf.full_file2).st_size == 0): # are there any variants in the final file? 
						VPOT_conf.update_existing_variants(VPOT_conf.working_file1,"1") # No - then create it based on sample file 
					elif (os.stat(VPOT_conf.working_file1).st_size == 0): # are there any variants for this sample? 
						copyfile(VPOT_conf.full_file2,VPOT_conf.working_file1) # copy the new full file
						VPOT_conf.update_existing_variants(VPOT_conf.working_file1,"0") # then for existing variants, just add 0 to the end for this sample 
					else : # yes, then incorporate them 
						VPOT_conf.incorporate_this_src_into_full_file() #
					copyfile(VPOT_conf.full_file2,VPOT_conf.full_file1) # copy the new full file
					
					
				else : #
					print ("Issue with input file 2: ", line) #
					return 1 # error

###########################################################################################################
#
###########################################################################################################
def setup_for_this_src_file_1(source_vcf,file_line): #
#
	for src_line in source_vcf: # work each line of source vcf file 
#			print (src_line) #
		src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
		if ("##" not in src_line1[0]): #
#				print src_line1 #
			if ("#CHROM" in src_line1[0]): # find sample location
				for i, content in enumerate(src_line1): # return the value and index number of each item in the line array 
#								print "content-",content,"/",i				#
					if (content == file_line[1]) : # when filtering for sample ID 
						VPOT_conf.sample_loc=i #save sample location
						VPOT_conf.sample_ID=file_line[1] # save sample ID
#							print "SAMPLE: ",sample_loc #
					elif (content == "INFO") : # 
						VPOT_conf.INFO_loc=i #save sample location
#							print "INFO_loc: ",INFO_loc #
					elif (content == "FORMAT") : # 
						VPOT_conf.FORMAT_loc=i #save sample location
#							print "FORMAT_loc: ",FORMAT_loc #
					#	
				if ("#CHROM" not in VPOT_conf.header_ln): # is this first time
#						print "setup header" #
					VPOT_conf.header_ln='\t'.join(src_line1[:VPOT_conf.INFO_loc+1])+tab+VPOT_conf.GENE_NAME # save heading up to INFO (using INFO_loc +1 to denote the stop field)
					VPOT_conf.blank_variant_ln=VPOT_conf.header_ln # save heading up to INFO (using INFO_loc +1 to denote the stop field)
#						VPOT_conf.header_ln='\t'.join(src_line1[:VPOT_conf.INFO_loc+1]) # save heading up to INFO (using INFO_loc +1 to denote the stop field)
#						VPOT_conf.blank_variant_ln='\t'.join(src_line1[:VPOT_conf.INFO_loc+1]) #
			else : # variants lines 
				FORMAT1=re.split(':',src_line1[VPOT_conf.FORMAT_loc]) # split into file location and sample id
#				print ("format : ",FORMAT1) #
				for i, content in enumerate(FORMAT1): # return the value and index number of each item in the line array 
#					print ("content-",content,"/",i)				#
					for j in range(len(VPOT_conf.sample_coverage)): #
#						print (VPOT_conf.sample_coverage[j],"/",content) #
						if (content == VPOT_conf.sample_coverage[j]) : 	# look for FORMAT field 
							VPOT_conf.sample_coverage_loc[j]=i 			# save sample location	
							break #
#					print ("coverage : ",VPOT_conf.sample_coverage_loc,i) #
				source_vcf.close() # finish with source vcf file 
				return 0 # have setup all location values - ok to go back
#							print INFO1 #
#							print FORMAT1 #
#							print SAMPLE1 #

#ed	source_vcf.close() # finish with source vcf file 
	return 1 # not a clean exit
#
###########################################################################################################
#
###########################################################################################################
def setup_for_this_src_file(file_line): #
#
#	global GT_loc #
	
#
#	print "setup_for_this_src_file): ", file_line #
#	print info_msg2 # parameter file supplied then use it
	VPOT_conf.sample_loc=-1 #
	VPOT_conf.sample_ID=""
	VPOT_conf.INFO_loc=-1 #
	VPOT_conf.FORMAT_loc=-1 #
#	VPOT_conf.sample_coverage_loc = [-1,-1,-1,-1,-1] # location of VCF format codes for sample 
#	VPOT_conf.sample_coverage_loc = [-1,-1,-1,-1,-1,-1] # location of VCF format codes for sample 
	#
#	with open(file_line[0],'r', encoding="utf-8") as source_vcf : #
##
	return_val=0 #
	if file_line[0].endswith('.gz'):
		with gzip.open(file_line[0],'rt') as source_vcf : #
			if ( setup_for_this_src_file_1(source_vcf,file_line) != 0 ): #
				return_val=1 #
	else:
		with open(file_line[0],'r',encoding="utf-8") as source_vcf : #
			if ( setup_for_this_src_file_1(source_vcf,file_line) != 0 ): #
				return_val=1 #
#
	return return_val #
#
###########################################################################################################
#
###########################################################################################################
def work_this_src_file_1(source_vcf, wrkf1): #
#
##
#	print "work_this_src_file(file_line): " #
#	print working_file1 #
		for src_line in source_vcf: # work each line of source vcf file 
			src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
#			if ("#" not in src_line1[0]): # skip the header lines
			if (("#" not in src_line1[0]) and ("_" not in src_line1[0])): # skip the header lines and any alternate chromosome alignments variants
#				print (src_line1) #
#				print src_line1[VPOT_conf.sample_loc] #
# variants lines 
				Sample_coverage=0 # initial sample coverage total 
				Sample_GQ=0 # initial sample genotype quality score 
				SAMPLE1=re.split(':',src_line1[VPOT_conf.sample_loc]) # split the sample's FORMAT fields 
#				print sample_coverage_loc[1],"/", sample_coverage_loc[2] #
#				print SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]],"/", SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]] #
# have a coverage check on the sample 
#				print ("MaxCOverage : ",VPOT_conf.Maxcoverage) #
#				print ("Hete_balance : ",VPOT_conf.Hete_Balance) #
#				print ("Genotype Quality : ",VPOT_conf.Genotype_Quality) #
#				print "coverage_loc : ",VPOT_conf.sample_coverage_loc #
#				print (SAMPLE1) 
#				print ("GT: ",VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val],":",SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]]) 
#				print ("NR: ",VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val],"/",SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] )
#				print ("NV: ",VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val],"/",SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]] )
#				print ("DP: ",VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val],"/",SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]]) 
# check to see if we want this variant line for this sample
				if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "./.") and (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "0/.") and (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "./0") and (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] != "0/0") : # a valid alternate genotype 
#					print ("DP: ",VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]) 
#    				check coverage depth
					if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] == ".") : # no NR_val 
						SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] = "0"  # set it as zero 
					if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]] == ".") : # no NV_val 
						SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]] = "0"  # set it as zero 
					if (VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val] != -1 ) : #this sample have a coverage depth
						if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]] == ".") : # no DP_val 
							SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]] = "0"  # set it as zero 
						Sample_coverage=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]]) # save DP value
						Alt_reads=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]])/2 # save DP value
#						print ("DP") #
					if (VPOT_conf.sample_coverage_loc[VPOT_conf.AD_val] != -1 ) : #this sample have an allele depth
						if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.AD_val]] != ".") : # has AD_val 
							AD_values=re.split(',',SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.AD_val]]) # get the alleles depth
							if (len(AD_values) > 1 ) : #  
								Alt_reads=int(AD_values[1]) # save alternate read count
							else : # this AD only has alternate count
								Alt_reads=int(AD_values[0]) # save alternate read count
#						print ("AD") #
					if (VPOT_conf.sample_coverage_loc[VPOT_conf.GQ_val] != -1 ) : #this sample have a genotype quality score
						if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GQ_val]] == ".") : # no GQ_val 
							SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GQ_val]] = "0"  # set it as zero 
						Sample_GQ=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GQ_val]]) # save GQ value
#						print ("GQ") #
					if (VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val] != -1 ) and (VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val] != -1 ) : #this sample have a coverage depth from NR and NV
						if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] == ".") : # no NR_val 
							SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]] = "0"  # set it as zero 
						if (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]] == ".") : # no NV_val 
							SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]] = "0"  # set it as zero 
#						Sample_coverage=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]])+int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]]) # save DP value
						Sample_coverage=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]]) # save DP value
						Alt_reads=int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NV_val]]) # save alt reads value
#						print ("NR+NV") #
#					print ("TOT: ",str(Sample_coverage)) # 
#					print ("ALT_READ: ",str(Alt_reads)) # 
#					print ("pass") #
#					if (((VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val] == -1) or 
#						((VPOT_conf.is_number(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]])) and (int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.NR_val]]) >= int(VPOT_conf.Maxcoverage)))) and  #  NR
#						((VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val] == -1) or 
#						((VPOT_conf.is_number(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]])) and (int(SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.DP_val]]) >= int(VPOT_conf.Maxcoverage))))) : # DP
					VPOT_conf.QC_PASS=False #
#	QC checks
#					if (Sample_coverage >= int(VPOT_conf.Maxcoverage)) : 
					if ((Sample_coverage >= int(VPOT_conf.Maxcoverage)) and (Sample_GQ >= int(VPOT_conf.Genotype_Quality))) : 
						if ( Sample_coverage == 0 ) : # can't have divide by 0
							Sample_coverage=1 # set a dummay value
						if (int((Alt_reads/Sample_coverage)*100) >= int(VPOT_conf.Hete_Balance)) : # Pas QC for coverage and balance
							VPOT_conf.QC_PASS=True # Yes
#							print ("QC_PASS",VPOT_conf.QC_PASS) #
					GT_values=re.split('/',SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]]) # get the genotype fields
#				print GT_values #
					for j in range(len(GT_values)) : #
#					print (GT_values[j]) #
						if ( GT_values[j] not in VPOT_conf.Non_alt_GT_types ) : # when filtering for QC value 
#							print ("keep this variant1") #
							check_this_variant(src_line, wrkf1) #
							break # get out of for loop (GT_values)
#				else : #  NR
#					print ("not adding") #
#
#	source_vcf.close() # finish with source vcf file 
###########################################################################################################
#
###########################################################################################################
#
def work_this_src_file(file_line): #
#
##
#	print "work_this_src_file(file_line): " #
#	print working_file1 #
	wrkf1=open(VPOT_conf.working_file1,'w', encoding="utf-8") # 
#
	if file_line[0].endswith('.gz'):
		with gzip.open(file_line[0],'rt') as source_vcf : #
			work_this_src_file_1(source_vcf,wrkf1) #
	else:
		with open(file_line[0],'r',encoding="utf-8") as source_vcf : #
			work_this_src_file_1(source_vcf,wrkf1) #
#
	wrkf1.close() # finish with the output file 
#
#	print "sort unique" #
	COMMAND="sort -V -u -k1,5 "+VPOT_conf.working_file1+" > "+VPOT_conf.sort_file1 #  
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
#	print (src_line1) #
#	print (SAMPLE1) #
	if (population_frequency(src_line1[VPOT_conf.INFO_loc]) == 0 ): # check if variant with in filter
#		print ("saving this variant1") #
		VPOT_conf.gene_ref="NONE" #
		find_gene_ref(src_line1[VPOT_conf.INFO_loc])    # check if variant found != means yes 
#
		if (not VPOT_conf.QC_PASS) : # sample did not pass QC  
			gtype="." #
		elif (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]] == "1/1") : # a homozygous alt genotype 
			gtype="2" #
		else : #
			gtype="1" 
#		print (SAMPLE1[VPOT_conf.sample_coverage_loc[VPOT_conf.GT_val]],"-",gtype)#
		outline='\t'.join(src_line1[:VPOT_conf.INFO_loc+1])+tab+VPOT_conf.gene_ref+tab+gtype+nl #
		wrkf1.write(outline) #
#	else : #
#		print "did not pass PF" #
#
###########################################################################################################
#
###########################################################################################################
def population_frequency(info_ln): #
#
#
#	print "population_frequency(info_ln): " #
	#
	val=0 #
	INFO1=re.split(';|=',info_ln) # split into file location and sample id
#	print (INFO1) #
#	print (VPOT_conf.PF_array) #
	for j in range(len(VPOT_conf.PF_array)): #
#		print (VPOT_conf.PF_array[j][1]) # move to pred_array slot
		for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#			print ("content-",content,"/",i)				#
#			print (VPOT_conf.PF_array[j][1]) # move to pred_array slot
			if (content == VPOT_conf.PF_array[j][1]) : # the population freq annotation we want? 
				if ( INFO1[i+1] != ".") :	 # if numeric check, else ok as most likely "." to state no annotation 
					temp_val=INFO1[i+1]
				else : # not number, so set a default as no annotation in population_frequency check would mean novel	
					temp_val=-999
#				print ("temp_val :",temp_val) #
#				if ( VPOT_conf.is_number(INFO1[i+1]) ) : # numeric
				#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
#					if ( float(INFO1[i+1]) > float(VPOT_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
				if ( VPOT_conf.is_number(temp_val) ) : # numeric
				#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
#					print ("temp_val yes :",temp_val) #
#					print (VPOT_conf.PF_array[j][3]) #
					if ( VPOT_conf.PF_array[j][3] == "LE" ) :	 # lower or each to PF limit 
						if ( float(temp_val) > float(VPOT_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
							val=1 # do not want this variant 
							break #
					else : # GT  look for greater or equal to PF limit	
						if ( float(temp_val) < float(VPOT_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
							val=1 # do not want this variant 
							break #
	#
#	print ("val :",val)
	return val #
#	
###########################################################################################################
#
###########################################################################################################
def find_gene_ref(info_ln): #
#
#
#	print "find_gene_ref(info_ln): " #
	#
	val=0 #
	INFO1=re.split(';|=',info_ln) # split into file location and sample id
#	print INFO1 #
#	print PF_array #
#	print PF_array[j][1] # move to pred_array slot
	for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#		print "content-",content,"/",i				#
		if (content == VPOT_conf.GN_value) : # the population freq annotation we want? 
			VPOT_conf.gene_ref=INFO1[i+1] # save gene name 
			break #
#
###########################################################################################################
#
###########################################################################################################
def score_the_variants(): #
##
#	print "score_the_variants(): " #
	#
	with open(VPOT_conf.full_file1,'r',encoding="utf-8") as variants_file, open(VPOT_conf.working_file1,'w',encoding="utf-8") as score_file : # 
		for line1 in variants_file: # work each line of new sample vcf file 
			priority_score=0 # initialise score 
			line_parts=re.split('\t|\n|\r',line1) # split the variant up
#			print "line part 0 : ",line_parts[0] #
			if ("#CHROM" != line_parts[0]): #
#				print src_line1 #
#				
				if (len(VPOT_conf.PD_array) > 0 ) : #  
					priority_score=prioritise_variants_by_predictors(line_parts[VPOT_conf.INFO_loc]) # get priority score
				#
				if (len(VPOT_conf.VT_array) > 0 ) : #   
#				print "check VT" #
					priority_score=priority_score+prioritise_variants_by_VT_types(line_parts[VPOT_conf.INFO_loc]) # 
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
def prioritise_variants_by_predictors(INFO_details): #
#
#	global PD_array 
#
#	print "prioritise_variants_by_predictors(INFO_details): " #
#
	type=2 #
	a_value=a_score=0 #
#
	val=0 #
	INFO1=re.split(';|=',INFO_details) # split into file location and sample id
#	print INFO1 #
	for j in range(len(VPOT_conf.PD_array)): #
#		print VT_array[j][1] # move to pred_array slot
		for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#			print "content-",content,"/",i				#
			if (content == VPOT_conf.PD_array[j][1]) : # the predictor annotation we want? 
				limit=len(VPOT_conf.PD_array[j]) # 
#				print "PD_array: ", limit, VPOT_conf.PD_array[j] #
#				print a_score, "/", a_value #
				a_value=type+1 #
				a_score=a_value+1 #
				if (VPOT_conf.PD_array[j][type] == "A") : # character field 
					while (a_value < len(VPOT_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
						if ( INFO1[i+1] == VPOT_conf.PD_array[j][a_value] ) :		# is this the variant type we are looking for   
							val=val+int(VPOT_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
							break 										# finish 
						a_value+=2	 # bump to next value option   
						a_score+=2	 #   
				else : # numeric field 
#					if ( INFO1[i+1] != ".") :	 # if annotation then continue   
					if ( VPOT_conf.is_number(INFO1[i+1]) ) : # numeric
#						print INFO1[i+1] #
						while (a_value < len(VPOT_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#							print VPOT_conf.PD_array[j][a_value] #
#							print "INFO1:", INFO1[i+1] #
							if ( float(INFO1[i+1]) < float(VPOT_conf.PD_array[j][a_value])) :	 # is this the variant type we are looking for   
								val=val+int(VPOT_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
								break 										# finish 
#							work1=a_value+2 # end of this PD limits  
#							print "work1:", work1 # end of this PD limits  
#							if (work1 >= len(VPOT_conf.PD_array[j])-1) : # end of this PD limits  
							if ( a_value+2 >= len(VPOT_conf.PD_array[j])-1) : # end of this PD limits  
#								print ("checking last value") #
								if ( float(INFO1[i+1]) > float(VPOT_conf.PD_array[j][a_value])) :	 # check the last value, which is a greater then check   
									val=val+int(VPOT_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
							a_value+=2	 # bump to next value option   
							a_score+=2	 # 
							
#						print "a_score:", a_score #
#						print "len:", len(VPOT_conf.PD_array[j]) #
				break # done for this predictor annotation 
#
	return val #
#	
##
#	print info_msg2 # parameter file supplied then use it
###########################################################################################################
#
###########################################################################################################
def prioritise_variants_by_VT_types(INFO_details): #
#
#	global VT_array 
#
#	print "prioritise_variants_by_VT_types(INFO_details): " #
#	print VT_array #
	val=0 #
	INFO1=re.split(';|=',INFO_details) # split into file location and sample id
#	print INFO1 #
	for j in range(len(VPOT_conf.VT_array)): #
#		print VT_array[j][1] # move to pred_array slot
		for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#			print "content-",content,"/",i				#
			if (content == VPOT_conf.VT_array[j][1]) : # the variant functional annotation we want? 
				if ( INFO1[i+1] == VPOT_conf.VT_array[j][2]) :	 # is this the variant type we are looking for   
					val=val+int(VPOT_conf.VT_array[j][3])	 # yes - return the score for the variant 
#
	return val #
##
#
