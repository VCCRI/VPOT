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
Sample_loc=-1 #
#
info_msg1_1="VPOT stats : Invalid number of inputs, must have two :" 
info_msg1_2="VPOT stats : 1) output destination directory + prefix" #
info_msg1_3="VPOT stats : 2) Input file - output from VPOT prioritisation process" 
info_msg1_4="VPOT stats : 3) percentile value - output gene breakdown for variant above this percentile (default is all variants)" 
#
Ast_ln=VPOT_conf.nl+"*********************************************************************************************"+VPOT_conf.nl #
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
#	print (sys.argv) #
	supplied_args=len(sys.argv) #
#	print (supplied_args) #
#
	if (supplied_args > 5 ):  # arg [0] is the python program
		print (info_msg1_1+VPOT_conf.nl+info_msg1_2+VPOT_conf.nl+info_msg1_3+VPOT_conf.nl+info_msg1_4) #
		return 1 #
	else :
		VPOT_conf.output_dir=sys.argv[2] #
		VPOT_conf.input_file=sys.argv[3] #
		VPOT_conf.parameter_file=50 # set default
#		print (supplied_args)
		if (supplied_args > 4 ):  # there is a percentile supplied
			if (VPOT_conf.is_number(sys.argv[4])): # numeric
				VPOT_conf.parameter_file=sys.argv[4] #
				if ((int(VPOT_conf.parameter_file) < 1) or (int(VPOT_conf.parameter_file) > 99)): # not a valid percentile
					VPOT_conf.parameter_file=50 # set default to all
#		print (VPOT_conf.parameter_file)
	#
		VPOT_conf.final_output_file=VPOT_conf.output_dir+"variant_statistic_file_"+suffix+".txt" #
		VPOT_conf.working_file1=VPOT_conf.output_dir+"working_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.working_file2=VPOT_conf.output_dir+"working_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.full_file1=VPOT_conf.output_dir+"full_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.full_file2=VPOT_conf.output_dir+"full_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.sort_file1=VPOT_conf.output_dir+"sort_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.sort_file2=VPOT_conf.output_dir+"sort_file2_"+suffix+"_tmp.txt" #
		print ("output : ",VPOT_conf.final_output_file) #
	
	return 0 #
#
#
###########################################################################################################
def Total_variants_stats(Output_file): #
#
	global scorefile #
	global gene_loc #
	global infile_shape #
#
	Outln=Ast_ln+VPOT_conf.nl+"TOTAL STATS FOR - "+VPOT_conf.input_file+VPOT_conf.nl #
#
	j=scorefile.shape # this is the array dimension for scoring variants
	num_var=infile_shape[0]-1 # number of variants in input file - assume a header
	#
	VPOT_conf.txt_start = 3 # samples id starts here
	VPOT_conf.txt_end = infile_shape[1] # -1 for the blank element 
	VPOT_conf.Sample_ids = scorefile[0,VPOT_conf.txt_start:VPOT_conf.txt_end] #
#	print (VPOT_conf.Sample_ids)
	num_samples = k = VPOT_conf.txt_end-VPOT_conf.txt_start #
#	print (j)
#
	Outln=Outln+"Total number of samples : "+str(num_samples)+VPOT_conf.nl #
#
#	gene_loc=10 # where gene name filed is located
	condition = np.unique( scorefile[:,gene_loc] ) #
	Outln=Outln+"Total number of genes : "+str(len(condition)-1)+VPOT_conf.nl #
	#	
	Outln=Outln+"Total number of variants : "+str(num_var)+VPOT_conf.nl #
#
	
	Outln=Outln+"Total number of non-scoring variants : "+str(infile_shape[0]-j[0])+VPOT_conf.nl #
	Outln=Outln+"Total number of scoring variants : "+str(j[0]-1)+VPOT_conf.nl+VPOT_conf.nl #
#
# determine breakdown of percentage of ranking with the scorefile (variants that have scores)
	rank=scorefile[:,0] # take the ranking column
	rank[0]=-1 # replace the header element
	rank=rank.astype(np.float) #
	#
	percentage=np.array([0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01]) #
	for i in range(len(percentage)): #
		condition=np.where(rank >= percentage[i]) #
		condition=condition[0]
		pc=len(condition)#
		pc_val=percentage[i]*100 #
		Outln=Outln+"Variants in the "+str(int(pc_val))+" percentile in ranking - "+str(pc)+VPOT_conf.nl #
#	
	print(Outln)
	Output_file.write(Outln) # write the line to final output file 
	#
###########################################################################################################
#
#
###########################################################################################################
def Samples_stats(Output_file,sample,pos): #
#
# input file for filtering is the output from the VPOT process. This means the variant gene name is in a column named GENE_NAME.
#
	global scorefile #
	global gene_loc #
#	
	Outln=Ast_ln+VPOT_conf.nl+"STATS for Sample "+sample+" : "+VPOT_conf.nl+VPOT_conf.nl #
	#
#
#	process to create a subset for sample that only contains variants that occurs in the sample
#
	condition = np.where( scorefile[:,pos] != "0" ) # look for variants that have is not 0 for sample
	condition = condition[0] #
	index = [0] # remove the header line
	con1 = np.delete(condition, index) #
#	print (con1)
	Outln=Outln+"Total number of variants : "+str(len(con1))+VPOT_conf.nl #
#
	n1 = con1.tolist() # convert the row numbers for the sample variants (where !="0")
	sample_all_arr = scorefile[np.ix_(n1,)] # subset based on these rows as index from scorefile
#
	unique_genes = np.unique( sample_all_arr[:,gene_loc] ) # we want to know the total number of genes for all the variants in the sample
#	print (unique_genes)
	Outln=Outln+"Total number of genes : "+str(len(unique_genes))+VPOT_conf.nl #
#	sample_arr = scorefile[np.ix_(n1,)] #
#
#	process to create a subset for sample that only contains variants that have ranking value > than the given parameter
#
	rank=sample_all_arr[:,0] # take the ranking column
	rank=rank.astype(np.float) # convert it into floating point value
	percentile=float(VPOT_conf.parameter_file)/100 # setup the percentile argument given as a floating point value
#	print(str(percentile)) #
#	print(rank)
	condition=np.where(rank >= percentile) # find the variants that meets the percentile criteria
	condition = condition[0] #
#	print(condition)
	n1 = condition.tolist() #convert the row values into a list
#	print(n1)
	sample_arr = sample_all_arr[np.ix_(n1,)] # now subset the sample's variants to one that meet the ranking percentile criteria
#	print(sample_arr)
#
# variants per gene for ones that meet the percentile criteria
#
	Outln=Outln+VPOT_conf.nl+"Variants within Gene breakdown based on supplied " #
	if (VPOT_conf.parameter_file==0 ):  # is this the default percentile
		Outln=Outln+"1 percentile parameter :"+VPOT_conf.nl # yes
	else :
		Outln=Outln+str(VPOT_conf.parameter_file)+" percentile parameter :"+VPOT_conf.nl # a supplied value
#
	unique_genes = np.unique( sample_arr[:,gene_loc] ) #
	for i in range(len(unique_genes)): #
		sample_gene = np.where( sample_arr[:,gene_loc] == unique_genes[i] ) #
		sample_gene = sample_gene[0] #
		Outln=Outln+VPOT_conf.nl+"Number of variants in "+unique_genes[i]+" : "+str(len(sample_gene)) #
		#		print (sample_gene) #
#
#	print(Outln) #
	Output_file.write(Outln) # write the line to final output file 
	#
###########################################################################################################
#
###########################################################################################################
	
def file_lines_numpy(working_fn): #
#
	global scorefile #
#	print "Place input into array : " #
#
	infile = np.genfromtxt(working_fn,comments=None,delimiter="\t",dtype=np.str) # this is the full file
	i=infile.shape #
#
	condition = np.where( infile[:,1] != "0" ) #
	condition = condition[0] #
	n1 = condition.tolist() #
	scorefile = infile[np.ix_(n1,)] # create an array with only scoring variants
#	print (VPOT_conf.Sample_ids) #
#	print (VPOT_conf.txt_start, VPOT_conf.txt_end)
#	print (infile[0,])
#	print(i)#
	return i #
##
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	global scorefile #
	global gene_loc #
	global infile_shape # contains the array dimension of the input file, so you can get number of variants and number of samples
#
	VPOT_conf.init() #
	#
	print ("Variant Statistics - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
	#
	COMMAND="cut --complement -f3-10 "+VPOT_conf.input_file+" > "+VPOT_conf.working_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
# Now filter the input file by gene list 
	infile_shape = file_lines_numpy(VPOT_conf.working_file1) # set up the numpy array
	gene_loc = 2 # standard col location for gene name
	with open(VPOT_conf.final_output_file,'w',encoding="utf-8") as Output_file : # 
		Total_variants_stats(Output_file) #
		i = VPOT_conf.txt_start # where the samples start
		j = 0 #
		while (i < VPOT_conf.txt_end): # for number of samples
#			print (i)
			Samples_stats(Output_file,VPOT_conf.Sample_ids[j],i) #
			i = i+1
			j = j+1
	#
#
#
############################################################################################################