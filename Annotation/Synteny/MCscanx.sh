#!/bin/bash

##########################################################
#run step wise
##########################################################

agat_convert_sp_gff2bed.pl --gff ficifolium_final_annotation_sorted.gff3 -o C._ficifolium.agatconvert.bed
agat_convert_sp_gff2bed.pl --gff cq_genomic.gff3 -o cq_genomic.bed
#These bed files were manually converted into the 4 column version that MCscanX wanted and dups removed - chr gene start stop
#this was done in excel by inserting column 4 into position 2
#Thus we went from C._ficifolium.agatconvert.gff3 to C._ficifolium.agatconvert.bed .
#This was cleaned, removing all -T1, -T2, etc and special characters, plus making everything lowercase.
#This created C._ficifolium.agatconvert.cleaned.bed
#I also renamed all CqA1 to qa1, etc. keeping with the lowercase two letter standard

#Later runs of this software used GFFToMcBedFormat.py to convert the fic gff3 into the format mcscan wants
python GFFToMcBedFormat.py > FicAnnotationMCConvert.bed

#reformat the quinoa fasta so headers are only gene names
bash ReformatFasta.sh

bash formatcdsprot.sh

blastp -query ficifolium.proteins.putative_function.fasta -subject 10.1.protein.cleaned.faa -outfmt 6 -evalue 1e-10 -max_hsps 5 -max_target_seqs 5 -out cfxbv.blast
#manually remove all -t1 -t2 -RA etc.
awk '{print $1, $2, $12}' cfxbv.blast > cfxbv.homology


mkdir CfBv
#move into dir
cd CfBv
#manually move and cat bed files together
cat *.bed > Cf_Bv.bed


tools/MCScanX/MCScanX_h Cf_Bv
