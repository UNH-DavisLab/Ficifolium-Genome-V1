#!/bin/bash

threads=22

#run the following line once for fai generation 
#samtools faidx V6CircularizedPolishedMitogenomenextpolish.fa
ref="V6CircularizedPolishedMitogenomenextpolish.fa"

freebayes-parallel <(fasta_generate_regions.py "$ref.fai" 200000) $threads \
-f $ref --gvcf --genotype-qualities --use-best-n-alleles 4 -g 5000 F2/PIA.MT.bam F2/QIA.MT.bam F2/3.MT.bam F2/4.MT.bam F2/5.MT.bam F2/6.MT.bam > RawMT.vcf

# F2/PIA.MT.bam F2/QIA.MT.bam F2/1.MT.bam F2/2.MT.bam F2/3.MT.bam F2/10-2.MT.bam F2/10.MT.bam F2/11-2.MT.bam F2/11.MT.bam F2/12-2.MT.bam F2/1-2.MT.bam F2/12.MT.bam F2/13-2.MT.bam F2/13.MT.bam F2/14-2.MT.bam F2/14.MT.bam F2/15-2.MT.bam F2/15.MT.bam F2/16-2.MT.bam F2/16.MT.bam F2/17-2.MT.bam F2/17.MT.bam F2/18-18.MT.bam F2/18-2.MT.bam F2/18.MT.bam F2/21.MT.bam F2/2-2.MT.bam F2/23.MT.bam F2/24.MT.bam F2/31.MT.bam F2/3-2.MT.bam F2/33.MT.bam F2/35.MT.bam F2/37.MT.bam F2/39.MT.bam F2/41.MT.bam F2/4-2.MT.bam F2/42.MT.bam F2/43B.MT.bam F2/44.MT.bam F2/4.MT.bam F2/5-2.MT.bam F2/5.MT.bam F2/6-2.MT.bam F2/6.MT.bam F2/7-2.MT.bam F2/7.MT.bam F2/8-2.MT.bam F2/8.MT.bam F2/9-2.MT.bam F2/9.MT.bam 
