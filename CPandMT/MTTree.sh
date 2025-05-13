#!/bin/bash

threads=78

#cat *.fa > allseqs.fasta

mafft --thread $threads --adjustdirectionaccurately AllSeqs.fasta > allseqs.msa.fasta
raxmlHPC -T $threads -s allseqs.msa.fasta -n AllSeqsRAX -m GTRCAT -p 64295 -o Beta_vulgaris-CQRefBased
