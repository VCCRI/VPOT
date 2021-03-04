# python3 convert_vep_multiple_transcript_vcf_to_official_vcf_format.py -i invcf.vcf -o outvcf.vcf
# developer - Emma Rath / VCCRI - 2021
#
# This program reads in a vcf that has been annotated with VEP, and may or may not have other INFO and may or may not have samples.

# An example input VEP VCF header line:
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|TSSDistance|miRNA|LoFtool|MPC|existing_InFrame_oORFs|existing_OutOfFrame_oORFs|existing_uORFs|five_prime_UTR_variant_annotation|five_prime_UTR_variant_consequence|LoF|LoF_filter|LoF_flags|LoF_info">

# An example input line (split for readability):
# chrX	1392611	.	C	T	668.60	PASS	
# AC=1;AF=0.5;AN=2;BaseQRankSum=1.1;ClippingRankSum=0;DP=59;EMQ=60;QD=11.94;ReadPosRankSum=-1.32;SOR=1.876;VQSLOD=20.69;culprit=MQRankSum;ANNOVAR_DATE=2018-04-16;
# Func.refGene=upstream;Gene.refGene=SLC25A6;GeneDetail.refGene=dist\x3d498;Func.wgEncodeGencodeBasicV33=upstream;Gene.wgEncodeGencodeBasicV33=SLC25A6;
# HGMDPRO2019p4_RANKSCORE=.;CLINVAR20201010_ALLELEID=.;CLINVAR20201010_CLNVI=.;BAYESDEL_ADDAF=.;BAYESDEL_NOAF=.;avsnp150=rs71254418;
# LBJ26ALL_SiPhy_29way_logOdds=.;regsnp_fpr=.;regsnp_disease=.;regsnp_splicing_site=.;tmcsnpdb=.;spliceai_site=.;spliceai_acceptor_gain_delta_score=.;
# VC_het_cnt=40;VC_hom_cnt=8;VC_freq=0.05957;VARIANT_TYPE=upstream_downstream;ALLELE_END;ANNOVAR_DATE=2018-04-16;CADD_RawScore=0.588078;CADD_PHRED=7.396;ALLELE_END;
# CSQ=T|upstream_gene_variant|MODIFIER|SLC25A6|ENSG00000169100|Transcript|ENST00000381401|protein_coding|||||||||||498|-1||HGNC|HGNC:10992||Ensembl|||||||498||0.263||||||||||,T|upstream_gene_variant|MODIFIER|LINC00106|ENSG00000236871|Transcript|ENST00000430235|lncRNA|||||||||||4414|1||HGNC|HGNC:31843||Ensembl|||||||4414||||||||||||,T|upstream_gene_variant|MODIFIER|LINC00106|ENSG00000236871|Transcript|ENST00000434938|lncRNA|||||||||||3816|1||HGNC|HGNC:31843||Ensembl|||||||3816||||||||||||,T|upstream_gene_variant|MODIFIER|SLC25A6|ENSG00000169100|Transcript|ENST00000475167|processed_transcript|||||||||||508|-1||HGNC|HGNC:10992||Ensembl|||||||508||0.263||||||||||,T|upstream_gene_variant|MODIFIER|SLC25A6|ENSG00000169100|Transcript|ENST00000484026|processed_transcript|||||||||||523|-1||HGNC|HGNC:10992||Ensembl|||||||523||0.263||||||||||,T|upstream_gene_variant|MODIFIER|SLC25A6|293|Transcript|NM_001636.4|protein_coding|||||||||||498|-1||EntrezGene|HGNC:10992||RefSeq|||||||498||0.263||||||||||,T|upstream_gene_variant|MODIFIER|LINC00106|751580|Transcript|NR_130733.1|lncRNA|||||||||||4414|1||EntrezGene|HGNC:31843||RefSeq|||||||4414||||||||||||,T|missense_variant|MODERATE|LOC105373102|105373102|Transcript|XM_011546186.1|protein_coding|1/3||||263|245|82|S/L|tCa/tTa|||1||EntrezGene|||RefSeq|||||||||||||||||||,T|non_coding_transcript_exon_variant|MODIFIER|LOC105373102|105373102|Transcript|XR_951284.1|misc_RNA|1/4||||263|||||||1||EntrezGene|||RefSeq|||||||||||||||||||,T|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR00000243759|promoter||||||||||||||||||||||||||||||||||||,T|TF_binding_site_variant|MODIFIER|||MotifFeature|ENSM00207005835|||||||||||||-1|||||||ENSPFM0006|8|N|0.025|HOXB2::TBX21&MGA::DLX2&MGA::DLX3&HOXB2::EOMES&HOXB2::TBX3&MGA::EVX1|||||||||||||	
# GT:AD:DP:GQ:HP:PL:PQ	0/1:31,25:56:99:1392561-2,1392561-1:676,0,870:1418.36

# An example output line (split for readability):
# chrX	1392611	.	C	T	668.60	PASS	
# AC=1;AF=0.5;AN=2;BaseQRankSum=1.1;ClippingRankSum=0;DP=59;EMQ=60;QD=11.94;ReadPosRankSum=-1.32;SOR=1.876;VQSLOD=20.69;culprit=MQRankSum;ANNOVAR_DATE=2018-04-16;
# Func.refGene=upstream;Gene.refGene=SLC25A6;GeneDetail.refGene=dist\x3d498;Func.wgEncodeGencodeBasicV33=upstream;Gene.wgEncodeGencodeBasicV33=SLC25A6;
# HGMDPRO2019p4_RANKSCORE=.;CLINVAR20201010_ALLELEID=.;CLINVAR20201010_CLNVI=.;BAYESDEL_ADDAF=.;BAYESDEL_NOAF=.;avsnp150=rs71254418;
# LBJ26ALL_SiPhy_29way_logOdds=.;regsnp_fpr=.;regsnp_disease=.;regsnp_splicing_site=.;tmcsnpdb=.;spliceai_site=.;spliceai_acceptor_gain_delta_score=.;
# VC_het_cnt=40;VC_hom_cnt=8;VC_freq=0.05957;VARIANT_TYPE=upstream_downstream;ALLELE_END;ANNOVAR_DATE=2018-04-16;CADD_RawScore=0.588078;CADD_PHRED=7.396;ALLELE_END;
# VEP_Allele=T|T|T|T|
# VEP_Consequence=upstream_gene_variant|upstream_gene_variant|missense_variant|TF_binding_site_variant|
# VEP_IMPACT=MODIFIER|MODIFIER|MODERATE|MODIFIER|
# VEP_SYMBOL=SLC25A6|LINC00106|LOC105373102|
# VEP_Gene=ENSG00000169100|ENSG00000236871|105373102||
# VEP_Feature_type=Transcript|Transcript|Transcript|MotifFeature|
# VEP_Feature=ENST00000381401|ENST00000430235|XM_011546186.1|ENSM00207005835|
# VEP_BIOTYPE=protein_coding|lncRNA|protein_coding||
# VEP_EXON=||1/3||
# VEP_INTRON=||||
# VEP_HGVSc=||||
# VEP_HGVSp=||||
# VEP_cDNA_position=||263||
# VEP_CDS_position=||245||
# VEP_Protein_position=||82||
# VEP_Amino_acids=||S/L||
# VEP_Codons=||tCa/tTa||
# VEP_Existing_variation=||||
# VEP_DISTANCE=498|4414|||
# VEP_STRAND=-1|1|1|-1|
# VEP_FLAGS=||||
# VEP_SYMBOL_SOURCE=|HGNC|HGNC|EntrezGene||
# VEP_HGNC_ID=HGNC:10992|HGNC:31843||
# VEP_REFSEQ_MATCH=||||
# VEP_SOURCE=Ensembl|Ensembl|||
# VEP_REFSEQ_OFFSET=||RefSeq||
# VEP_MOTIF_NAME=|||ENSPFM0006|
# VEP_MOTIF_POS=|||8|
# VEP_HIGH_INF_POS=|||N|
# VEP_MOTIF_SCORE_CHANGE=|||0.025|
# VEP_TRANSCRIPTION_FACTORS=|||HOXB2::TBX21&MGA::DLX2&MGA::DLX3&HOXB2::EOMES&HOXB2::TBX3&MGA::EVX1|
# VEP_TSSDistance=498|4414||
# VEP_miRNA=||
# VEP_LoFtool=0.263||||
# VEP_MPC=||||
# VEP_existing_InFrame_oORFs=||||
# VEP_existing_OutOfFrame_oORFs=||||
# VEP_existing_uORFs=||||
# VEP_five_prime_UTR_variant_annotation=||||
# VEP_five_prime_UTR_variant_consequence=||||
# VEP_LoF=||||
# VEP_LoF_filter=||||
# VEP_LoF_flags=||||
# VEP_LoF_info=||||

import sys
import csv
import subprocess
import argparse
import VPOT_conf #


def parse_arguments(args):
    parser = argparse.ArgumentParser(description='Reformat VEP INFO in VCF file to conform to VCF specification. Each VEP field will become an INFO field. Values of multiple transcripts will be separated by bar in each INFO field.')
    parser.add_argument('-i', dest='invcf', 
                        help='Input VCF annotated with VEP')
    parser.add_argument('-o', dest='outvcf', 
                        help='Output VCF file within reformatted VEP INFO fields')
    return parser.parse_args()


def read_invcf(input_file_path):
  with open(input_file_path) as invcf_handle:
    inlines = invcf_handle.readlines()
  inlines_without_trailing_carriage_return = []
  for inline in inlines:
    new_inline = inline[ 0 : len(inline)-1 ]
    inlines_without_trailing_carriage_return.append(new_inline)
  return inlines_without_trailing_carriage_return


def extract_and_convert_and_write_out_header(invcf, outvcf):

  # input line:
  # ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|TSSDistance|miRNA|LoFtool|MPC|existing_InFrame_oORFs|existing_OutOfFrame_oORFs|existing_uORFs|five_prime_UTR_variant_annotation|five_prime_UTR_variant_consequence|LoF|LoF_filter|LoF_flags|LoF_info">

  # output lines:
  # ##INFO=<ID=VEP_Allele,Number=1,Type=String,Description="VEP Allele">
  # ##INFO=<ID=VEP_Consequence,Number=1,Type=String,Description="VEP Consequence">

  list_of_vep_hdrs = []
  id_of_old_vep_hdr = '##INFO=<ID=CSQ,'

  for inline in invcf:
    if (len(inline) >=1):
      char1 = inline[0:1]
      if (char1 == '#'):
        this_inline_id = ''
        if (len(inline) >= len(id_of_old_vep_hdr)):
          this_inline_id = inline[0:len(id_of_old_vep_hdr)]
        if (this_inline_id == id_of_old_vep_hdr):
          split1 = inline.split('"')
          extract1 = split1[1]
          split2 = extract1.split('Format: ')
          extract2 = split2[1]
          split3 = extract2.split('|')
          for this_vep_field in split3:
            this_vep_hdr = "VEP_" + str(this_vep_field)
            list_of_vep_hdrs.append(this_vep_hdr)
            outline = '##INFO=<ID=' + str(this_vep_hdr) + ',Number=1,Type=String,Description="VEP ' + str(this_vep_field) + '">' + "\n"
            outvcf.write(outline)
        else:
            outline = str(inline) + "\n"
            outvcf.write(outline)
  return list_of_vep_hdrs


def reformat_VEP_info_entry(old_info_field, list_of_vep_hdrs):

  # input format:
  # CSQ=T|upstream_gene_variant|MODIFIER|SLC25A6,T|missense_variant|MODERATE|LOC105373102,T|TF_binding_site_variant|MODIFIER|,T|missense_variant|MODERATE|LOC105373102

  # output format:
  # VEP_Allele=T|T|T|T;VEP_Consequence=upstream_gene_variant|upstream_gene_variant|missense_variant|TF_binding_site_variant;VEP_IMPACT=MODIFIER|MODIFIER|MODERATE||MODIFIER;VEP_SYMBOL=SLC25A6|LINC00106|LOC105373102|

  list_of_vep_info_fields = [''] * len(list_of_vep_hdrs)

  transcripts = old_info_field.split(',')
  for i in range( 0, len(transcripts) ): # process each transcript
    this_transcript = transcripts[i]
    fields = this_transcript.split('|')
    for j in range( 0, len(fields) ): # process the data fields in each transcript
      this_field = fields[j]
      if (i == 0): # processing first transcript
        this_info_field_output = str(list_of_vep_hdrs[j]) + '=' + str(this_field)
      else: # processing subsequent transcripts
        this_info_field_output = str(list_of_vep_info_fields[j]) + '|' + str(this_field)
      list_of_vep_info_fields[j] = this_info_field_output

  new_info_fields = ''
  for j in range( 0, len(list_of_vep_hdrs) ):
    this_info_field_output = list_of_vep_info_fields[j]
    bits = this_info_field_output.split('=')
    if (bits[1] == ''):
      this_info_field_output = this_info_field_output + '.'
    if (new_info_fields == ''):
      new_info_fields = this_info_field_output
    else:
      new_info_fields = new_info_fields + ';' + this_info_field_output

  return new_info_fields


def reformat_vcf_info(old_info, list_of_vep_hdrs):

  new_info = ''
  infields = old_info.split(";")
  for infield in infields:
    bits = infield.split('=')
    this_info_key = bits[0]
    if (this_info_key == 'CSQ'):
      this_info_data = bits[1]
      new_infield = reformat_VEP_info_entry(this_info_data, list_of_vep_hdrs)
    else:
      new_infield = str(infield)

    if (new_info == ''):
      new_info = new_infield
    else:
      new_info = new_info + ";" + new_infield
    if (new_info == ''):
      new_info = '.'
      
  return new_info


def reformat_vcf_records(invcf, outvcf, list_of_vep_hdrs):

  for inline in invcf:
    if (len(inline) >=1):
      char1 = inline[0:1]
      if (char1 != '#'):

        # process this input vcf data record
        infields = inline.split("\t")
        old_info = infields[7]

        # convert vep in info field from vep format to correct vcf format
        new_info = reformat_vcf_info(old_info, list_of_vep_hdrs)

        infields[7] = new_info

        # output the vcf record with correctly formatted info field
        outline = str(infields[0])
        for i in range(1, len(infields)):
          outline = outline + "\t" + str(infields[i])
        outline = str(outline) + "\n"
        outvcf.write(outline)
  return


def convert_vep_multiple_transcript_vcf_to_official_vcf_format(input_file_path, output_file_path):

  invcf = read_invcf(input_file_path)

  outvcf = open(output_file_path, 'w')

  list_of_vep_hdrs = extract_and_convert_and_write_out_header(invcf, outvcf)

  reformat_vcf_records(invcf, outvcf, list_of_vep_hdrs)

  return

#
###########################################################################################################
#
###########################################################################################################
#
def main():

##  args = parse_arguments(args)
  input_file_path = VPOT_conf.input_file
  output_file_path = VPOT_conf.final_output_file

  convert_vep_multiple_transcript_vcf_to_official_vcf_format(input_file_path, output_file_path)

  #print(' ')
  #print('End of convert_vep_multiple_transcript_vcf_to_official_vcf_format')

