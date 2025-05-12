#!/bin/bash

python3 quartet.py TeloExplorer -i ordered_new_ficifolium_corrected_PCsInLinearOrder.fasta -c plant 
python3 quartet.py CentroMiner -t 78 -p OrderedCorrected -i ordered_new_ficifolium_corrected_PCsInLinearOrder.fasta
