###########################################################################################################
# 
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import VPOT_conf #
from shutil import copyfile #
#
tab='\t' # 
nl='\n' #
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
info_msg1_1="VPOT merge : Invalid number of inputs, must have three :" 
info_msg1_2="VPOT merge : 1) output destination directory + prefix" #
info_msg1_3="VPOT merge : 2) Input file 1 - file to be merged" 
info_msg1_4="VPOT merge : 3) Input file 2 - file to be merged" 
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
	if (supplied_args != 5 ):  # arg [0] is the python program
		print (info_msg1_1+VPOT_conf.nl+info_msg1_2+VPOT_conf.nl+info_msg1_3+VPOT_conf.nl+info_msg1_4) #
		return 1 #
	else :
		VPOT_conf.output_dir=sys.argv[2] #
		VPOT_conf.input_file=sys.argv[3] #
		VPOT_conf.parameter_file=sys.argv[4] #
#		print (VPOT_conf.input_file)
#		print (VPOT_conf.parameter_file)
	#
		VPOT_conf.final_output_file=VPOT_conf.output_dir+"variant_merge_file_"+suffix+".txt" #
		VPOT_conf.working_file1=VPOT_conf.output_dir+"working_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.working_file2=VPOT_conf.output_dir+"working_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.full_file1=VPOT_conf.output_dir+"full_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.full_file2=VPOT_conf.output_dir+"full_file2_"+suffix+"_tmp.txt" #
		VPOT_conf.sort_file1=VPOT_conf.output_dir+"sort_file1_"+suffix+"_tmp.txt" #
		VPOT_conf.sort_file2=VPOT_conf.output_dir+"sort_file2_"+suffix+"_tmp.txt" #
		print ("output : ",VPOT_conf.final_output_file) #
	
	return 0 #
#
############################################################################################################
###########################################################################################################
#
###########################################################################################################
def merge_the_input(infile):
#
	#
#	print ("merge_the_input(): "+infile) #
	MEOF = True #
#
	try:
		os.remove(VPOT_conf.sort_file1) #
	except OSError:
		pass #
	
	with open(infile,'r',encoding="utf-8") as input_file, open(VPOT_conf.sort_file2,'r',encoding="utf-8") as master_file, open(VPOT_conf.sort_file1,'w',encoding="utf-8") as out_file : # 
		inl1 = input_file.readline().rstrip()
		inl1p = inl1.split('\t')
		header_line = inl1p
		h_cnt=len(header_line) #
#		print(str(header_line)+":"+str(h_cnt)) #
		emplst =["0"] * h_cnt
#		print(str(emplst)) #
#		print(str(inl1p[2:10])) #
		mstl1 = master_file.readline().rstrip()
		mstl1p = mstl1.split('\t')
		m_cnt=len(mstl1p) #
#		print(str(h_cnt)+"/"+str(m_cnt)) #
#		print(str(mstl1p[2:10])) #
#
		while MEOF : #do while we still have record in master
			if (inl1p[2:6] == mstl1p[2:6]):
#				print("match")
				out_file.write(mstl1+tab+inl1+nl) # write the line to final output file 
				inl1 = input_file.readline().rstrip()
				inl1p = inl1.split('\t')
#				print(str(inl1p[2:10])) #
				mstl1 = master_file.readline().rstrip()
				mstl1p = mstl1.split('\t')
#				print(str(mstl1p[2:10])) #
			else :	# then input line is not in master so need to pad master
#				print("not match")
				out_file.write(mstl1+tab+tab.join(emplst)+nl) # write the line to final output file 
				mstl1 = master_file.readline().rstrip()
				mstl1p = mstl1.split('\t')
			#
			if (mstl1 == ''):
				MEOF = False 
#
	c1=m_cnt+1 
	c2=c1+10 
	COMMAND="cut -f "+str(c1)+"-"+str(c2)+" --complement "+VPOT_conf.sort_file1+" > "+VPOT_conf.sort_file2 # create a working input file - by reducing to only the cols we need 
#	print (COMMAND) #
	subprocess.call(COMMAND, shell=True) # do it in shell
#
	return 0 #
#
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	global scorefile #
	global gene_loc #
	global infile_shape #
#
	VPOT_conf.init() #
	#
	print ("Merge Variant Output files - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
#
	#print(VPOT_conf.input_file+" "+VPOT_conf.parameter_file+" "+VPOT_conf.final_output_file+" "+VPOT_conf.working_file1)
	# create a working input file - by reducing to only the cols we need 

	COMMAND="tail -n+2 "+VPOT_conf.input_file+" > "+VPOT_conf.working_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="cut -f 1-11 "+VPOT_conf.working_file1+" > "+VPOT_conf.sort_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="tail -n+2 "+VPOT_conf.parameter_file+" > "+VPOT_conf.working_file2 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="cut -f 1-11 "+VPOT_conf.working_file2+" >> "+VPOT_conf.sort_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="head -n1 "+VPOT_conf.input_file+" | cut -f 1-11 > "+VPOT_conf.sort_file2 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="sort -u -k3,3V -k4,4n -k6,7 "+VPOT_conf.sort_file1+" >> "+VPOT_conf.sort_file2 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="head -n1 "+VPOT_conf.input_file+" > "+VPOT_conf.full_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="sort -u -k3,3V -k4,4n -k6,7 "+VPOT_conf.working_file1+" >> "+VPOT_conf.full_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="head -n1 "+VPOT_conf.parameter_file+" > "+VPOT_conf.full_file2 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="sort -u -k3,3V -k4,4n -k6,7 "+VPOT_conf.working_file2+" >> "+VPOT_conf.full_file2 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
	print ("merge_the_input(): "+VPOT_conf.input_file) #
	merge_the_input(VPOT_conf.full_file1) # merge input file 1 into the new 
	print ("merge_the_input(): "+VPOT_conf.parameter_file) #
	merge_the_input(VPOT_conf.full_file2) # merge input file 2 into the new 
#
	COMMAND="head -n1 "+VPOT_conf.sort_file2+ " > "+VPOT_conf.final_output_file # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="tail -n+2 "+VPOT_conf.sort_file2+" > "+VPOT_conf.working_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
	COMMAND="sort -u -k2,2nr -k3,3V -k4,4n -k6,7 "+VPOT_conf.working_file1+" >> "+VPOT_conf.final_output_file # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
############################################################################################################