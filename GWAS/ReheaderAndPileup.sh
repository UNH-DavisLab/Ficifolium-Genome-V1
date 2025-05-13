#!/bin/bash
#reheader all vcf files according to their directory. This will also call mpileup from bam files and produce a vcf for stacks pop pipe in the stacksv1linkage dir.

threads=78
ref="ficifolium/V1/ReferenceFiles/ordered_new_ficifolium_with_contigs_corrected.fa" #changed 3/16/23
reffai="ficifolium/V1/ReferenceFiles/ordered_new_ficifolium_with_contigs_corrected.fa.fai"


for dir in ficifolium/V1/F2toV1/*/
do
    Label=$(basename $dir)
    BamF=$(find $dir -name *.bam)

    FileNameBam="${Label}.bam"
    VCFFile="${Label}.g.vcf.gz"
    VarBam="${Label}.variant.bam"
    echo "${dir}${Label}.bam"
    gatk AddOrReplaceReadGroups -I "${dir}${Label}.bam" -O "${dir}${Label}.rehead.bam" -RGID "ReGr${Label}" -RGLB "lib1-${Label}" -RGPL ILLUMINA -RGPU unit1 -RGSM "${Label}"
    samtools index -b -@ $threads "${dir}${Label}.rehead.bam"
 
done

freebayes-parallel <(fasta_generate_regions.py $reffai 100000) $threads \
-f $ref --gvcf --genotype-qualities --use-best-n-alleles 4 -g 5000 1/1.rehead.bam 2/2.rehead.bam 4/4.rehead.bam 5/5.rehead.bam 6/6.rehead.bam 7/7.rehead.bam 8/8.rehead.bam 9/9.rehead.bam 10/10.rehead.bam 11/11.rehead.bam 12/12.rehead.bam 13/13.rehead.bam 14/14.rehead.bam 15/15.rehead.bam 16/16.rehead.bam 17/17.rehead.bam 18/18.rehead.bam 1-2/1-2.rehead.bam 2-2/2-2.rehead.bam 3-2/3-2.rehead.bam 4-2/4-2.rehead.bam 5-2/5-2.rehead.bam 6-2/6-2.rehead.bam 7-2/7-2.rehead.bam 8-2/8-2.rehead.bam 9-2/9-2.rehead.bam 10-2/10-2.rehead.bam 11-2/11-2.rehead.bam 12-2/12-2.rehead.bam 13-2/13-2.rehead.bam 14-2/14-2.rehead.bam 15-2/15-2.rehead.bam 16-2/16-2.rehead.bam 17-2/17-2.rehead.bam 18-2/18-2.rehead.bam 18-18/18-18.rehead.bam 21/21.rehead.bam 23/23.rehead.bam 24/24.rehead.bam 31/31.rehead.bam 33/33.rehead.bam 35/35.rehead.bam 37/37.rehead.bam 39/39.rehead.bam 41/41.rehead.bam 42/42.rehead.bam 43B/43B.rehead.bam 44/44.rehead.bam | \
vcffilter -f "QUAL > 990" > F2sWithParentsfreebayes.vcf


#bgzip F2sWithParentsfreebayes.vcf > F2sWithParentsfreebayes.vcf.gz

