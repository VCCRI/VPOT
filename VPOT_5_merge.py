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
info_msg1_1="VPOT merge : Invalid number of inputs, must have two :" 
info_msg1_2="VPOT merge : 1) output destination directory + prefix" #
info_msg1_3="VPOT merge : 2) Input file containing list of VPOL files to be merged" 
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
	if (supplied_args != 4 ):  # arg [0] is the python program
		print (info_msg1_1+VPOT_conf.nl+info_msg1_2+VPOT_conf.nl+info_msg1_3) #
		return 1 #
	else :
		VPOT_conf.output_dir=sys.argv[2] #
		VPOT_conf.input_file=sys.argv[3] # file containing list of VPOLs
#		print (VPOT_conf.input_file)
	#
		VPOT_conf.final_output_file=VPOT_conf.output_dir+"variant_merge_file_"+suffix+".txt" #
		VPOT_conf.temp_output_file=VPOT_conf.output_dir+"variant_merge_file_"+suffix+"_tmp.txt" #
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
def setup_for_merge(file1):
#
#
#	print ("adding input file to full master : ",file1) #
	COMMAND="tail -n+2 "+file1+" > "+VPOT_conf.working_file1 # create a working input file for file1- 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="cut -f 1-11 "+VPOT_conf.working_file1+" >> "+VPOT_conf.sort_file1 # reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
#	COMMAND="sort -u -k3,3V -k4,4n -k5,7 "+VPOT_conf.sort_file1+" > "+VPOT_conf.sort_file2 # create a working input file - by reducing to only the cols we need 
	COMMAND="sort -u -k3,3V -k4,4n -k6,7 "+VPOT_conf.sort_file1+" > "+VPOT_conf.sort_file2 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
	copyfile(VPOT_conf.sort_file2,VPOT_conf.sort_file1) # 
	
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
	
	COMMAND="head -n1 "+infile+" > "+VPOT_conf.full_file1 # create a working input file 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="tail -n+2 "+infile+" > "+VPOT_conf.working_file1 # create a working input file for file1- 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="sort -u -k3,3V -k4,4n -k6,7 "+VPOT_conf.working_file1+" >> "+VPOT_conf.full_file1 # create a working input file 
	subprocess.call(COMMAND, shell=True) # do it in shell
	
	with open(VPOT_conf.full_file1,'r',encoding="utf-8") as input_file, open(VPOT_conf.temp_output_file,'r',encoding="utf-8") as master_file, open(VPOT_conf.sort_file1,'w',encoding="utf-8") as out_file : # 
		inl1 = input_file.readline().rstrip()
		inl1p = inl1.split('\t')
		header_line = inl1p
		h_cnt=len(header_line) #
#		print(str(header_line)+":"+str(h_cnt)) #
		emplst =["0"] * h_cnt
#		print(str(emplst)) #
#		print(str(inl1p[2:6])) #
		mstl1 = master_file.readline().rstrip()
		mstl1p = mstl1.split('\t')
		m_cnt=len(mstl1p) #
#		print(str(mstl1p)+":"+str(m_cnt)) #
#		print(str(h_cnt)+"/"+str(m_cnt)) #
#		print(str(mstl1p[2:6])) #
#
		while MEOF : #do while we still have record in master
#			print("comp-",str(inl1p[2:7])) #
#			print("mst-",str(mstl1p[2:7])) #
#
#			if (inl1p[2:7] == mstl1p[2:7]):
			if ((inl1p[2:4] == mstl1p[2:4]) and (inl1p[5:7] == mstl1p[5:7])):
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
	COMMAND="cut -f "+str(c1)+"-"+str(c2)+" --complement "+VPOT_conf.sort_file1+" > "+VPOT_conf.temp_output_file # create a working input file - by reducing to only the cols we need 
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
	with open(VPOT_conf.input_file,'r',encoding="utf-8") as check_file : # 
		inVPOLs = alist = [line.rstrip() for line in check_file]
	print(inVPOLs) #
	if (len(inVPOLs) < 2 ): #
		print("No need to merge as only one input file is provided.") #
	else : #
		#
# setup the master file (VPOT_conf.sort_file2)
		for x in range(0, len(inVPOLs)): #
			setup_for_merge(inVPOLs[x]) # merge input file into the new file (VPOT_conf.sort_file2)
			#
		COMMAND="head -n1 "+inVPOLs[0]+" | cut -f 1-11 > "+VPOT_conf.temp_output_file #  
		subprocess.call(COMMAND, shell=True) # do it in shell
		COMMAND="cat "+VPOT_conf.sort_file2+" >> "+VPOT_conf.temp_output_file # complete the creation of the merge master file  
		subprocess.call(COMMAND, shell=True) # do it in shell
		#
#		print ("merge_the_input(): "+inVPOLs[0]) #
		for x in range(0, len(inVPOLs)): #
			print ("merge_the_input(): "+inVPOLs[x]) #
			merge_the_input(inVPOLs[x]) # merge input file into the new file (VPOT_conf.sort_file2)
#
	COMMAND="head -n1 "+VPOT_conf.temp_output_file+ " > "+VPOT_conf.final_output_file # setup the final output file with the header 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="tail -n+2 "+VPOT_conf.temp_output_file+" > "+VPOT_conf.working_file1 # take the variants portion 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
	COMMAND="sort -u -k2,2nr -k3,3V -k4,4n -k6,7 "+VPOT_conf.working_file1+" >> "+VPOT_conf.final_output_file # sort the output based on ranking score 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
############################################################################################################