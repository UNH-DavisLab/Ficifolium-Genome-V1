#!/bin/bash
#wrapper for vcf2gwas
#conda activate vcf2gwas
#################################################################################
#NOTE: THE VCF FILE IN THIS DIRECTORY HAS BEEN MODIFIED BY THE FOLLOWING LINES 14-24
#Once the vcf has been modified, do not run any lines above vcf2gwas
#################################################################################
VCFFile="F2sWithParentsfreebayes.vcf"
PhenoFile="OrigPheno.csv"
Ref="ficifolium/V1/ReferenceFiles/ordered_new_ficifolium_with_contigs_corrected.fa"
GFF="ficifolium/V1/funannotate/FinalFiles/ficifolium_masked_final_annotation.maker.integrated.gff.csv"
Threads=22

#First we need to swap the infoline "NUMALT" for something else
#sed -i 's/NUMALT/NUMALTERNATE/g' $VCFFile
#Next we normalize the vcf by splitting all multiallelic sites into biallelic SNPs
#bcftools norm --multiallelics -any --threads $Threads -f $Ref $VCFFile -o temp.vcf
#mv temp.vcf $VCFFile

#rename the individuals in the file because having the "-" in the name field of the vcf breaks it. "-" replaced with "r". This file will likely need to be edited if there are further revisions.
#bcftools reheader --samples VCFRename.txt $VCFFile -o temp.vcf
#mv temp.vcf $VCFFile

#bgzip -k $VCFFile 

#-nq for no quality checking. 

vcf2gwas -v "$VCFFile.gz" -nq -pf $PhenoFile -gf $GFF -p 1 -lmm 4 -o First21/FlowTime
vcf2gwas -v "$VCFFile.gz" -nq -pf $PhenoFile -gf $GFF -p 2 -lmm 4 -o First21/Height
vcf2gwas -v "$VCFFile.gz" -nq -pf $PhenoFile -gf $GFF -p 3 -lmm 4 -o First21/NumBranch
vcf2gwas -v "$VCFFile.gz" -nq -pf $PhenoFile -gf $GFF -p 4 -lmm 4 -o First21/BranchAngle
vcf2gwas -v "$VCFFile.gz" -nq -pf $PhenoFile -gf $GFF -p 5 -lmm 4 -o First21/IntLen
