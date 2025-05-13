#!/bin/bash

threads=22

#run the following line once for fai generation 
ref="CpPath2.fasta"
ls="LSC.fasta"
ss="SSC.fasta"

#samtools faidx $ref
#samtools faidx $ls
#samtools faidx $ss

#F2 numbers 1, 2, and 3 correspond to samples 4, 5, and 6, with 3 being the F1. PIA is the female parent, QIA the male. 

freebayes-parallel <(fasta_generate_regions.py "$ref.fai" 200000) $threads \
-f $ref --gvcf --genotype-qualities F2/PIA.CP.bam F2/QIA.CP.bam F2/3.CP.bam F2/4.CP.bam F2/5.CP.bam F2/6.CP.bam > RawCP.vcf

#All possible samples below:
# F2/PIA.CP.bam F2/QIA.CP.bam F2/1.CP.bam F2/2.CP.bam F2/3.CP.bam F2/10-2.CP.bam F2/10.CP.bam F2/11-2.CP.bam F2/11.CP.bam F2/12-2.CP.bam F2/1-2.CP.bam F2/12.CP.bam F2/13-2.CP.bam F2/13.CP.bam F2/14-2.CP.bam F2/14.CP.bam F2/15-2.CP.bam F2/15.CP.bam F2/16-2.CP.bam F2/16.CP.bam F2/17-2.CP.bam F2/17.CP.bam F2/18-18.CP.bam F2/18-2.CP.bam F2/18.CP.bam F2/21.CP.bam F2/2-2.CP.bam F2/23.CP.bam F2/24.CP.bam F2/31.CP.bam F2/3-2.CP.bam F2/33.CP.bam F2/35.CP.bam F2/37.CP.bam F2/39.CP.bam F2/41.CP.bam F2/4-2.CP.bam F2/42.CP.bam F2/43B.CP.bam F2/44.CP.bam F2/4.CP.bam F2/5-2.CP.bam F2/5.CP.bam F2/6-2.CP.bam F2/6.CP.bam F2/7-2.CP.bam F2/7.CP.bam F2/8-2.CP.bam F2/8.CP.bam F2/9-2.CP.bam F2/9.CP.bam 
