#!/bin/bash

#align and make tree from cp sequences
#manually move beta vulgaris to the top of the file to set it as an outgroup by default in iqtree

#conda activate iqtree

mafft --thread 22 allseqs.fasta > allseqs.msa.fasta
raxmlHPC -T 22 -s allseqs.msa.fasta -n AllSeqsRAX -m GTRCAT -o OU343016.1_Beta_vulgaris
