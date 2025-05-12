#!/bin/bash

##########################################################################################################
# SETUP AND CONFIGURATION
##########################################################################################################
# Activate the required environment
mamba activate funannotate

Threads=23
Genome="ordered_new_ficifolium_corrected.fa.softmasked.fasta"
IsoseqFasta="ficifolium_flnc_polished.hq.fasta"
OrthoDB="Viridiplantae.fa"

##########################################################################################################
# PREPARE ENVIRONMENT AND ALIGN ISOSEQ TO GENOME
##########################################################################################################

# Build Braker3 container if necessary
singularity build braker3_lr.sif docker://teambraker/braker3:isoseq

# Align IsoSeq reads to genome and convert SAM to BAM
minimap2 -t${Threads} -ax splice:hq -uf $Genome $IsoseqFasta > isoseq.sam
samtools view -bS --threads ${Threads} isoseq.sam -o isoseq.bam

##########################################################################################################
# RUN BRAKER FOR GENE PREDICTION
##########################################################################################################

# Copy Augustus config to fix new species write error
singularity exec -B $PWD:$PWD braker3_lr.sif cp -R /opt/Augustus/config .

# Run Braker for annotation using IsoSeq and protein homology
singularity exec -B ${PWD}:${PWD},config:/augustus_config braker3_lr.sif \
braker.pl --gff3 --AUGUSTUS_CONFIG_PATH=/augustus_config --genome=$Genome --prot_seq=$OrthoDB \
--bam=isoseq.bam --threads=${Threads} --species=C.ficifolium --busco_lineage=embryophyta_odb10

##########################################################################################################
# QUALITY CHECK WITH BUSCO
##########################################################################################################

# Run BUSCO on the predicted proteome and coding sequences
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.3.2_cv1 \
busco -i braker/braker.aa -o ProteomeBrakerEmbryophyta_odb10 -l embryophyta_odb10/ -m proteins -c ${Threads}

docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.3.2_cv1 \
busco -i braker/braker.codingseq -o TranscriptomeBrakerEmbryophyta_odb10 -l embryophyta_odb10/ -m transcriptome -c ${Threads}

# Both produced 98.5% completion

##########################################################################################################
# GET BASE ANNOTATION STATS
##########################################################################################################

# Get annotation stats on BRAKER annotation
agat_sp_statistics.pl --gff braker/braker.gff3 > BrakerStats.txt

##########################################################################################################
# RUN MAKER FOR AED SCORES AND REFINEMENT
##########################################################################################################

maker

# Merge GFFs from different annotations
agat_sp_merge_annotations.pl -o joinedMaker.gff3 \
-f C2/CA/Cf1/Cf1.gff \
-f 4A/59/Cf2/Cf2.gff \
-f 2F/39/Cf3/Cf3.gff \
-f 1F/90/Cf4/Cf4.gff \
-f 3A/20/Cf5/Cf5.gff \
-f B9/83/Cf6/Cf6.gff \
-f 10/DA/Cf7/Cf7.gff \
-f F4/6A/Cf8/Cf8.gff \
-f 65/5E/Cf9/Cf9.gff

# Remove fasta and rename annotation file - manual step
mv joinedMakerNoFasta.gff3 ficifolium.gff3

cat C2/CA/Cf1/Cf1.maker.proteins.fasta \
4A/59/Cf2/Cf2.maker.proteins.fasta \
2F/39/Cf3/Cf3.maker.proteins.fasta \
1F/90/Cf4/Cf4.maker.proteins.fasta \
3A/20/Cf5/Cf5.maker.proteins.fasta \
B9/83/Cf6/Cf6.maker.proteins.fasta \
10/DA/Cf7/Cf7.maker.proteins.fasta \
F4/6A/Cf8/Cf8.maker.proteins.fasta \
65/5E/Cf9/Cf9.maker.proteins.fasta > maker.proteins.fasta

cat C2/CA/Cf1/Cf1.maker.transcripts.fasta \
4A/59/Cf2/Cf2.maker.transcripts.fasta \
2F/39/Cf3/Cf3.maker.transcripts.fasta \
1F/90/Cf4/Cf4.maker.transcripts.fasta \
3A/20/Cf5/Cf5.maker.transcripts.fasta \
B9/83/Cf6/Cf6.maker.transcripts.fasta \
10/DA/Cf7/Cf7.maker.transcripts.fasta \
F4/6A/Cf8/Cf8.maker.transcripts.fasta \
65/5E/Cf9/Cf9.maker.transcripts.fasta > maker.transcripts.fasta

##########################################################################################################
# QUALITY CHECKS
##########################################################################################################

# Run BUSCO on the predicted proteome and coding sequences
docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.3.2_cv1 \
busco -i maker.proteins.fasta -o MakerProteomeEmbryophyta_odb10 -l embryophyta_odb10/ -m proteins -c 22

docker run -u $(id -u) -v $(pwd):/busco_wd ezlabgva/busco:v5.3.2_cv1 \
busco -i maker.transcripts.fasta -o MakerTranscriptomeEmbryophyta_odb10 -l embryophyta_odb10/ -m transcriptome -c 22

# Both buscos produced 98.4% completeness

# Get annotation stats on final
agat_sp_statistics.pl --gff joinedMaker.gff3 > MakerStats.txt

# Find AED value spread
perl AED_cdf_generator.pl -b 0.025 joinedMakerNoFasta.gff3 > MakerAEDs.txt

##########################################################################################################
# FUNCTIONAL ANNOTATION: SWISS-PROT AND INTERPROSCAN
##########################################################################################################

# BLAST against UniProt Swiss-Prot
makeblastdb -in uniprot_sprot.fasta -blastdb_version 5 -title "uniprot_sprot" -dbtype prot
blastp -query maker.proteins.fasta -db uniprot_sprot.fasta -num_threads 72 -evalue 1e-6 \
-max_hsps 1 -max_target_seqs 1 -outfmt 6 -out ficifolium.proteins.sprot.blastp

# InterProScan for domain and GO term annotation
interproscan.sh -appl pfam -dp -cpu 22 -f TSV -goterms -iprlookup -pa -t p -i maker.proteins.fasta \
-o ficifolium.proteins.iprscan


##########################################################################################################
# RENAME IDS AND INTEGRATE GFF3
##########################################################################################################

# Rename IDs in GFF and FASTA files
maker_map_ids --prefix Cfic_ --justify 8 ficifolium.gff3 > ficifolium.map

# Run custom script to filter genes removed by maker from braker. Manually changed gene names to "Removed"
python UpdateMap.py

cp ficifolium.gff3 ficifolium.renamed.gff3
cp maker.proteins.fasta maker.proteins.renamed.fasta
cp maker.transcripts.fasta maker.transcripts.renamed.fasta
cp ficifolium.proteins.iprscan ficifolium.proteins.renamed.iprscan
cp ficifolium.proteins.sprot.blastp ficifolium.proteins.sprot.renamed.blastp

# Map renamed IDs to all annotation files
map_gff_ids ficifolium_complete_deleted.map ficifolium.renamed.gff3
map_fasta_ids ficifolium_complete_deleted.map maker.proteins.renamed.fasta 
map_fasta_ids ficifolium_complete_deleted.map maker.transcripts.renamed.fasta
map_data_ids ficifolium_complete_deleted.map ficifolium.proteins.renamed.iprscan
map_data_ids ficifolium_complete_deleted.map ficifolium.proteins.sprot.renamed.blastp

# Run custom script to update the old "g" gene names left within blast lines. Run it twice because chunking seems to exclude some entries.
python UpdateBlastGenes.py
mv ficifolium.renamedBlast.gff3 ficifolium.renamed.gff3
python UpdateBlastGenes.py
mv ficifolium.renamedBlast.gff3 ficifolium.renamed.gff3

# Functional annotation in FASTA files
maker_functional_fasta uniprot_sprot.fasta ficifolium.proteins.sprot.renamed.blastp maker.proteins.renamed.fasta > ficifolium.proteins.putative_function.fasta
maker_functional_fasta uniprot_sprot.fasta ficifolium.proteins.sprot.renamed.blastp maker.transcripts.renamed.fasta > ficifolium.cds-transcripts.putative_function.fasta

# Update GFF with functional annotations
maker_functional_gff uniprot_sprot.fasta ficifolium.proteins.sprot.renamed.blastp ficifolium.renamed.gff3 > ficifolium.codons_functional.gff
ipr_update_gff ficifolium.codons_functional.gff ficifolium.proteins.renamed.iprscan > ficifolium_masked_codons_final_annotation.integrated.gff3

##########################################################################################################
# CLEAN AND SORT
##########################################################################################################

# Clean
agat_convert_sp_gxf2gxf.pl -g ficifolium_masked_codons_final_annotation.integrated.gff3 -o ficifolium_final_annotated.cleaned.gff3

# Sort the final GFF - may not be required
gff3_sort -g ficifolium_final_annotated.cleaned.gff3 -og ficifolium_final_annotation_sorted.gff3

