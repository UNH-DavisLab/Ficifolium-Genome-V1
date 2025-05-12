#!/bin/bash

tools/RepeatModeler-2.0.2a/BuildDatabase -name V1Corrected -engine ncbi V1/ReferenceFiles/ordered_new_ficifolium_corrected.fa
tools/RepeatModeler-2.0.2a/RepeatModeler -database V1Corrected -pa 23 -LTRStruct -engine ncbi >& DenovoRepeats.out &

#Genome Masking for annotation
RepeatMasker -pa 79 -gff -a -html -lib consensi.fa.classified -dir MaskerOutput V1/ReferenceFiles/ordered_new_ficifolium_corrected.fa

#Files moved manaully to correct directory.
tools/RepeatMasker/util/./calcDivergenceFromAlign.pl -s ordered_new_ficifolium_corrected.fa.align.divsum ordered_new_ficifolium_corrected.fa.align
#genome size for v1 is 729997426 
tools/RepeatMasker/util/./createRepeatLandscape.pl -g 729997426 -div ordered_new_ficifolium_corrected.fa.align.divsum > ordered_new_ficifolium_corrected.html
