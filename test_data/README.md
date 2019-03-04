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

 - ExAC
 - gnomAD
 - 1000G (2015aug)
 - Polyphen2_HDIV
 - PolyPhen2_HVAR
 - CADD
 - MCAP
 - SIFT
 - LRT
 - Mutation Assessor
 - Mutation Taster
 - MetaSVM
 - GERP++
 - PhyloP20way_mammalian
 - SiPhy_29way
      
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