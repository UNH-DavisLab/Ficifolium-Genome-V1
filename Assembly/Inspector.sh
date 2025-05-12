#!/bin/bash
#inspector pipeline for genome cleaning
#https://github.com/Maggi-Chen/Inspector
#conda activate inspector
tools/Inspector/inspector.py -t 75 -c V1/ReferenceFiles/ordered_new_ficifolium.fa -r V1/hifireads/m54336U_220422_175202.hifi_reads.fastq.gz -o Ordered_new_fic/ --datatype hifi 
tools/Inspector/inspector-correct.py -t 75 -i Ordered_new_fic/ --datatype pacbio-hifi 
