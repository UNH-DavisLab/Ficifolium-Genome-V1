#!/bin/bash
#this script aligns these F2 reads to the fic mito and cp genome, then moves them to their
#respective locations in another directory

threads=22
CPBaseDir="ficifolium/V1/cytoplasm/CP/F2"
MTBaseDir="ficifolium/V1/cytoplasm/MT/F2"
#The following ref files are the most up-to-date as of 1-18-24, embplant path 2 and mito v6.
CPRef="ficifolium/V1/cytoplasm/CP/CpPath2.fasta"
MTRef="ficifolium/V1/cytoplasm/MT/V6CircularizedPolishedMitogenomenextpolish.fa"

bwa index $CPRef
bwa index $MTRef

#alignment
for dir in /home/clayton/bioinformatics/ficifolium/V1/F2toV1/*/
do
	FileName=$(basename $dir)
    read1=$(find $FileName -name "*R1*")
    read2=$(find $FileName -name "*R2*")
    echo $read1
    echo $read2

    #Alignment and sorting
    echo "BWA aligning CP Using ${read1} and ${read2}"
    bwa mem -t $threads $CPRef $read1 $read2 > "$FileName.CP.sam"	#align reads to ref
    samtools view -@ $threads -b -F 4 "$FileName.CP.sam" > "$FileName.CP.bam"	#convert sam to bam
    samtools sort -@ $threads "$FileName.CP.bam" > "$FileName.CP.sorted.bam"	#sort bam
    gatk AddOrReplaceReadGroups -I "$FileName.CP.sorted.bam" -O "$FileName.CP.rehead.bam" -RGID "ReGr${FileName}" -RGLB "lib1-${FileName}" -RGPL ILLUMINA -RGPU unit1 -RGSM "${FileName}"
    mv "$FileName.CP.rehead.bam" "$FileName.CP.bam"
    rm "$FileName.CP.sorted.bam"
    rm "$FileName.CP.sam"
    samtools index -b -@ $threads "$FileName.CP.bam"
    samtools consensus -@ $threads -a "$FileName.CP.bam" -o "$FileName.CP.fa"


    echo "BWA aligning MT Using ${read1} and ${read2}"
    bwa mem -t $threads $MTRef $read1 $read2 > "$FileName.MT.sam"       #align reads to ref
    samtools view -@ $threads -b -F 4 "$FileName.MT.sam" > "$FileName.MT.bam"   #convert sam to bam
    samtools sort -@ $threads "$FileName.MT.bam" > "$FileName.MT.sorted.bam"    #sort bam
    gatk AddOrReplaceReadGroups -I "$FileName.MT.sorted.bam" -O "$FileName.MT.rehead.bam" -RGID "ReGr${FileName}" -RGLB "lib1-${FileName}" -RGPL ILLUMINA -RGPU unit1 -RGSM "${FileName}"
    mv "$FileName.MT.rehead.bam" "$FileName.MT.bam"
    rm "$FileName.MT.sorted.bam"
    rm "$FileName.MT.sam"
    samtools index -b -@ $threads "$FileName.MT.bam"
    samtools consensus -@ $threads -a "$FileName.MT.bam" -o "$FileName.MT.fa"



    #cleanup and move files to correct directory

    cp "$FileName.CP.bam" "$CPBaseDir/$FileName.CP.bam"
    mv "$FileName.CP.bam" "${dir}/$FileName.CP.bam"
    cp "$FileName.CP.fa" "$CPBaseDir/$FileName.CP.fa"
    mv "$FileName.CP.fa" "${dir}/$FileName.CP.fa"
    cp "$FileName.CP.bam.bai" "$CPBaseDir/$FileName.CP.bam.bai"
    mv "$FileName.CP.bam.bai" "${dir}/$FileName.CP.bam.bai"

    cp "$FileName.MT.bam" "$MTBaseDir/$FileName.MT.bam"
    mv "$FileName.MT.bam" "${dir}/$FileName.MT.bam"
    cp "$FileName.MT.fa" "$MTBaseDir/$FileName.MT.fa"
    mv "$FileName.MT.fa" "${dir}/$FileName.MT.fa"
    cp "$FileName.MT.bam.bai" "$MTBaseDir/$FileName.MT.bam.bai"
    mv "$FileName.MT.bam.bai" "${dir}/$FileName.MT.bam.bai"

done



#put together a VCF within the CP and MT directories
