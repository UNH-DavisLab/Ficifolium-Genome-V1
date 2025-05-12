#!/bin/bash

#Assemble
hifiasm -o ficifolium.asm -t 22 m54336U_220422_175202.hifi_reads.fastq.gz 

#Convert gfa to fasta
awk '/^S/{print ">"$2"\n"$3}' ficifolium.asm.bp.p_ctg.gfa | fold > ficifolium.asm.bp.p_ctg.fa
