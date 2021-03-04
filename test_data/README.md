# Files in VPOT test data directory

| File name  |  Description  | 
|:--:|:--:|
| FAM_1.vcf | Test family VCF ( Proband : PATIENT1, Father : PATIENT2, Mother : PATIENT3)  |
| FAM_1.txt | Test Family TXT ( Proband : PATIENT1, Father : PATIENT2, Mother : PATIENT3) |
| default_0.001_variants_parameters_PPF.txt | Default Proritisation Parameter File (PPF) |
| test_VCF_sample_list.txt | test input samples file - contain a list of input VCFs and corresponding sample IDs  |
| test_TXT_sample_list.txt | test input samples file - contain a list of input TXTs and corresponding sample IDs
| test_gene.txt | test gene list - for use with genef option |
| test_sample_set.ped | test sample filter ped file - where PATIENT1 is affected, and PATIENT2 and PATIENT3 are unaffected - for use with samplef option |
| merge_t1.txt, merge_t2.txt, merge_t3.txt | test VPOL files - for use with merge option |
| test_merge_VPOL_list.txt | test input VPOL list file - contains a list of input VPOL files for merge option |

# Default Proritisation Parameter File (PPF) 
## 1. How was it derived?
 VPOT is designed to allow the users to customise their prioritisation process based on their disease of study and the annotations they have chosen; however, we have provided a default PPF (with a list of recommended annotations) based on our experience with the complex disease, Congenital Heart Disease.  The weighting criteria for scoring within the default PPF are set up to identify pathogenic prediction by in-silico predictors (e.g. CADD, Polyphen-2), with recommendation from literature, and to highlight the most disruptive variants such as stop-gain and splicing variants.
     
## 2. Required Minimum Annotations for use of default PPF

In order to use the default PPF it is required that these databases be included in your list of annotation databases. They cover in general most of the commonly used population frequecy databases as well as *in-silico* predictors. All of these annotation databases can be found in the ANNOVAR databases download.

| Annotation Database | Type  |
|--|--|
|ExAC  | Population Frequency  |
| gnomAD | Population Frequency |
| 1000G (2015aug) |  Population Frequency  |
| Polyphen2_HDIV | Functional Prediction |
| PolyPhen2_HVAR | Functional Prediction |
| CADD | Ensemble Score |
| MCAP | Ensemble Score |
| SIFT | Functional Prediction |
| LRT | Functional Prediction |
| Mutation Assessor | Functional Prediction |
| Mutation Taster | Functional Prediction |
| MetaSVM | Ensemble Score |
| GERP++ | Conservation Score |
| PhyloP20way_mammalian | Conservation Score |
| SiPhy_29way | Conservation Score |
 
## 3. Dealing with missing annotation/pathogenicity predictors values

In our CHD study we have post-annotated pathogenicity scores that lacks a value with a -999 value in the input VCFs. We then score these -999 values like any other pathogenicity score, as a middle ground between non-pathogenic and pathogenic.  This allows us to still highlight the variant above benign variants but do not eclipse true deleterious variants. 
     
# Tutorial on the use of VPOT 

## 1. Priority Option
**priority**  - prioritisation tool
this option performs the variant proritisation process on the input samples VCF files. It will score each variant found for the supplied samples based on the weighting affixed to the annotations supplied in the Prioritisation Parameter File, PPF.
     
  **Command line** : python3 VPOT.py priority <location for output file+prefix> < file of input VCF files> < parameter file> 

 Test example  :     
 - change to test_data directory
  - python3 ../VPOT.py priority testout_ test_VCF_sample_list.txt default_0.0001_variants_parameters_PPF.txt  

 result :    
 - testout_VPOL_final_output_file_XXXXXXXXXX.txt  
     
 
## 2. Gene Filtering  Option   
**genef**     - gene filter
   this option performs variant filtering of the VPOL based on genes supplied as input.

   **Command line** : python3 VPOT.py genef <location for output file+prefix> < VPOT prioritiy output> < gene list>    
  
  Test example  :     
 - change to test_data directory
 - python3 ../VPOT.py samplef testout_ testout_final_output_file_XXXXXXXXXX.txt test_gene.txt   
 
 result :    
 - testout_gene_filtered_output_file_XXXXXXXXXX.txt  
     

## 3. Samples and Inheritance Model Filtering Option   
**samplef**   - samples and inheritance model filtering
   this option performs variant filtering of the VPOL based on a suplied ped format file. It can be used for a simple case-control filtering or an inheritance model filtering for a family trio.
  
**Inheritance Models** - 

 1.  De novo (**DN**) model - filters for variants that only exist in the proband and not in any of the parents. 
 2. Autosomal dominant (**AD**) model - filters for variants that exist in both the proband and affected parent but not in the unaffected parent. 
 3. Autosomal recessive (**AR**) model  - filter for variants that is a homozygous alternate in the proband and heterozygous in both parents. 
 4. Compound heterozygous (**CH**) model - provides a filter that returns variants in genes that have both proband-parental and proband-maternal specific variants.


 **For case-control samples filtering :** 
 **Command line** : python3 VPOT.py samplef < location for output file+prefix > < VPOT prioritiy output > < sample selection file >    

 Test example  :     
   - change to test_data directory
   - python3 ../VPOT.py samplef testout_samplef_ testout_final_output_file_XXXXXXXXXX.txt test_sample_set.ped   

 result  :  
  - testout_samplef_variant_filtered_output_file_XXXXXXXX.txt  
  
 **For inheritance model filtering :** 
 **Command line** : python3 VPOT.py samplef < location for output file+prefix > < VPOT prioritiy output > < sample selection file > < proband sample ID > < inheritance model - DN/AD/AR/CH >   

 Test example  :     
   - change to test_data directory
   - python3 ../VPOT.py samplef testout_samplef_ testout_final_output_file_XXXXXXXXXX.txt test_sample_set.ped  PATIENT1 AR 

 result  :  
  - testout_samplef_variant_filtered_output_file_PATIENT1_AR_XXXXXXXX.txt  

## 4. Statistic Report Option    
 **stats**     - general statistics on the VPOT priority output file
this option returns a summary statistic file for the VPOL supplied. It provide a small report listing the number of variants (the total number of variants, the number of scoring variants, the number of non-scoring variants), the number of genes and the number of samples. There is a breakdown of the number of variants found in each score 10% percentile range. The top 20 variants are also listed. A breakdown for each sample is also provided, with a table containing the number of variants in genes found above the percentile value supplied.

   **Command line** : python3 VPOT.py stats <location for output file+prefix> <VPOT prioritiy output> <a percentile value [1-99] for quick summary statistic>    
 
  for the input percentile value the summary provide a gene breakdown of variants that meet the percentile criteria for each sample
- to find the highest priority variant/gene use a high value eg 90 for 90th percentile   
- to list all variants/genes for the sample use 1 for 1st percentile    
- if no specific input value is given then a default value of 50 is used    

Test example  :     
   - change to test_data directory
 - python3 ../VPOT.py stats testout_ merge_t3.txt 50   

 result  :  
  - an output file contain a statistic summary report of samples/genes/variants found in VPOL merge_t3.txt 

## 5. Merge VPOLs Option    
**merge**     - merge multiple VPOT priority output list files into a single VPOT priority output list file
 
this option provide the ability to merge various VPOLs into one single VPOL. This function allows large cohorts to be split into small groups to speed up proritisation processing and then output to be re-consolidated back into one single large cohort VPOL for downstream analysis or filtering.     

   **Command line** : python3 VPOT.py merge <location for output file+prefix> < test input VPOL list file>    

Test example  :     
   - change to test_data directory
 - python3 ../VPOT.py merge testout_ merge_ test_merge_VPOL_list.txt   

 result  :  
  - a VPOL output file will be created containing all the variants and samples from the VPOLs inside the test_merge_VPOL_list.txt.  

## 6. Utility Option    
**utility**     - general VPOT utilties
 
this option provide the ability to use different utilities available in the VPOT program suite.

 | Utility name | Description  | command string |
|--|--|--|
|convertVEP  | convert a VEP annotated VCF to standard format VCF  | python3 VPOT.py utility convertVEP < location of input VEP annotated VCF file> < location of output standard VCF file>|

More details regarding each utility can be found in the input files setup section for "utility" option below.


# INPUT FILES SETUP

##  1. SAMPLES INPUT FILE


 FORMAT : (one sample per line):    
 
 location of VCF file/text file< tab >sample id   
 
**see test_VCF_sample_list.txt in the test_data directory for format example.**

##  2. VPOT PRIORITISATION PARAMETER FILE (PPF) FOR OPTION PRIORITY     

### 2.1. Setting up PF population filter in PPF

 To provide population frequency threshold for selection of variants.

 **FORMAT :**    PF	< population frequency annotation name in VCF>	< filter Value>

 Example :   PF	ExAC_ALL	0.01

 This will tell VPOT to use the ExAC_ALL annotation values as a variant filtering criteria. VPOT will return variants that are <= to the value given, in this case 0.01.
 Multiple PF lines can be provided if you want to filter based on a combination of population frequency datasets. Note it is a AND logical approach, so the return variant would have
 met all the PF criteria. 

**see default_0.001_variants_parameters_PPF.txt in the test_data directory for format example.**  

### 2.2. Setting up PD predictors in parameter file 

 For each predictor a range of predictor scores are given along with the corresponding VPOT value. This allows the user to determine the VPOT value to be assigned to a specific score or prediction.
 VPOT allows for as many breakdown/levels of differentiation as the user wants or the predictor needs.

 **FORMAT :**   
 
 PD	< predictor annotation name in VCF>	< (A)lpha/(N)umeric prediction type>	[ prediction value/score]	[VPOT value for prediction value/score] .....[repeat as many values as you need] 
 
 see example below or in parameter files in the test data folder. 

 Multiple PD lines can be provided if you want to use multiple predictors 

 **ALPHA PREDICTORS** 

 For predictors that return alphanumeric prediction categories the PD line will list each category and its assigned VPOT value.
 
 Example :
 MutationTaster predicts an alteration as one of four possible types:

 * disease causing - i.e. probably deleterious - D 
 * disease causing automatic - i.e. known to be deleterious, see section dbSNP / TGP / ClinVar / HGMD for details- A
 * polymorphism - i.e. probably harmless - P 
 * polymorphism automatic - i.e. known to be harmless, see section dbSNP / TGP / ClinVar / HGMD for details N

 PD	MutationTaster_pred	A	N	0	P	0	D	1	A	2	
 
 | Predictor Value | VPOT Value |
 |:---------------:|:----------:|
  N | 0
  P | 0
  D | 1 
  A | 2

 **NUMERIC PREDICTORS** 

 For predictors that return a numeric score the PD line will list scoring threshold and its assigned VPOT value.
 
 Example :
 CADD predicts an alteration with a score range from 0-30+ :

 * CADD < 10 - probably harmless 
 * CADD between 10 and 20 - median value for all possible canonical splice site changes and non-synonymous variants, as stated by CADD 
 * CADD > 20 has been recommend by some papers as a good pathogenicity threshold
 
 Example :   PD	CADD_phred	N	10	0	20	1	20	2
 
 | Predictor Value | VPOT Value |
 |:---------------:|:----------:|
  x < 10 | 0
  10 < x < 20 | 1
  20 < x | 2 
  
  so the parameter line can contain
 S1|V1|S2|V2|S3|V3|.......|Sn|Vn|Sn+1|Vn+1
 
 where S = the predictor score
 
  V = VPOT value
  
 Sn = Sn+1   
 
 | Predictor Value | VPOT Value |
 |:---------------:|:----------:|
  x < S1 | V1
  S2 < x < S3 | V2
  S3 < x < S4 | V3
  Sn < x < Sn+1 | Vn
  Sn+1 < x | Vn+1 

 **see default_0.001_variants_parameters_PPF.txt in the test_data directory for format example.**  

### 2.3. Setting up VT annotation in parameter file 

 To provide a point value to certain variant types that might not be well covered by predictors, eg STOPGAIN/SPLICING. 
 This allow highlighting of vertain variant types.

 **FORMAT :**    VT	< annotation field name in VCF that holds the variant type>	< VPOT_Value>

 Example :   VT	VARIANT_TYPE	exonic_stopgain_	50

 For variant type "exonic_stopgain" an extra 50 will be added to the variant's cumulative value.
 Multiple VT lines can be provided to provide different stratification of variants 

**see default_0.001_variants_parameters_PPF.txt in the test_data directory for format example.**  


### 2.4. Setting up GN gene symbol in parameter file 


 To provide the annotation field that contains the gene the variant is located in. 

 **FORMAT :**    GN	< annotation field name in VCF that holds the gene name>

 Example :   GN	Gene Symbol

 ONLY one field should be provided for this parameter option.

**see default_0.001_variants_parameters_PPF.txt in the test_data directory for format example.**  

### 2.5. Setting up QC threshold values in parameter file 


 To provide the Quality Control threshold for sample's genotype. 

 **FORMAT :**    QC	< Coverage | Hete_Balance > < threhold value>

 Example :   QC  Coverage   8

 For Coverage it is minumum number of reads that is acceptable.
 For Hete_Balance it is the minumum percentage acceptable between alternate 
 allele reads against total number of reads.
 
 A failed QC sample will be marked with "."
 
 Note: A variant will not be removed from the VPOL even if all samples fail the genotype QC threshold.

**see default_0.001_variants_parameters_PPF.txt in the test_data directory for format example.**  

### 2.6. Setting up VS threshold value in parameter file 

 To provide the Variant Score Threshold for the variants to output into the VPOL. 

 **FORMAT :**    VS	 Score < threhold value>

 Example :   VS Score   14

 This threshold value determine whether a variant will be included in the VPOL. It can be used to reduced the number of variants returned if the user want to set a minimum score. A variant is included in the VPOL if its Variant score is >= to this threshold.

**see default_0.001_variants_parameters_PPF.txt in the test_data directory for format example.**  

##  3. GENE SELECTION FILE FOR OPTION GENEF   
    
 The gene selection is based on a text file with a single gene name each line.
  
 **FORMAT :**  (one gene per line):     
 
 eg: ACTC1  

**see test_gene.txt in the test_data directory for format example.** 
   
##  4. SAMPLE SELECTION PED FILE FOR OPTION SAMPLEF   
    
 The sample selection is based on the pedigree ped file format, where the affected column is used to determine the selection of a variant.     
    
 If you are selecting variants based on one sample, eg all variants for a sample from a multi-sample VCF, then you can just specify that single sample in the 
 sample slection file, eg only want all sample PATIENT1 variants :
 
 |FID  |SAMPLE_ID |PAT |MAT |SEX |PHENOTYPE|
 |:---:|:--------:|:-------:|:-------:|:--:|:-------:|
 FAM_1|PATIENT1	|PATIENT2|PATIENT3|1|2  
    
 If you want variants that appear in PATIENT1 but not PATIENT2, a segregation request, then you will need a line for each sample 
 
 |FID  |SAMPLE_ID |PAT |MAT |SEX |PHENOTYPE|
 |:---:|:--------:|:-------:|:-------:|:--:|:-------:|
 FAM_1|PATIENT1	|PATIENT2|PATIENT3|1|2
 FAM_1|PATIENT2|ND1|ND2|1|1
   
 -  PATIENT1  2  = affected (variant IS in sample PATIENT1)
 
 -  PATIENT2  1  = unaffected (variant IS NOT in sample PATIENT2)
 
 so a combination of these values will determine if a variant is maintained or not. 
 For the above case a variant is maintain if it is found in PATIENT1 and not in PATIENT2. 
 
 This filtering option allows you to work with any number of samples whether they are related or not.
 
 Note: if there are more samples than the ones stated in the initial input file then they do not influence the variant selection. 

**see test_sample_set.ped in the test_data directory for format example.**

 ##  5. INPUT FILE FOR OPTION MERGE  
    
 **FORMAT :** (one VPOL per line): require a minimum of two files    
 
 location of VPOL file 
 
 **see test_merge_VPOL_list.txt in the test_data directory for format example.**

   ##  6. INPUT PARAMETERS FOR OPTION UTILITIES
   ### 6.1. convertVEP  

 This utility convert a VEP annotated VCF, which contain INFO annotation that are in a specific VEP format for multiple transcripts into a standard format VCF.  
 
   **Command line** : python3 VPOT.py utility convertVEP <location of input VEP annotated VCF file> < location of output standard VCF file>    

Test example  :     
   - change to test_data directory
 - python3 ../VPOT.py utility convertVEP merge test_vep_original.vcf test_vep_converted.vcf   

 result  :  
  - a VCF that follows the VCF format standard, which can be used as input to VPOT.  
 