#!/bin/bash
#this is designed to take the ficifolium assembly and align forward and reverse reads in each folder to it
#dont forget to index files used as ref. 

threads=50
RefLoc="ficifolium/V1/ReferenceFiles/ordered_new_ficifolium_with_contigs_corrected.fa" #changed 3/21/23

#to make the required ref files - the following two lines only need to be run once. 
#samtools faidx $RefLoc
#gatk CreateSequenceDictionary -R $RefLoc
#bwa index $RefLoc


#cat the fasta files together and remove old ones if there are duplicate R1 and R2s for a given sample

#for dir in ficifolium/V1/F2toV1/*/
#do
#    Label=$(basename $dir)
#    echo $dir
#    cat $dir/*R1* >> "${Label}Combined_R1.fastq.gz"
#    cat $dir/*R2* >> "${Label}Combined_R2.fastq.gz"
#    rm $dir/*001.fastq.gz
#    mv "${Label}Combined_R1.fastq.gz" $dir/"${Label}Combined_R1.fastq.gz"
#    mv "${Label}Combined_R2.fastq.gz" $dir/"${Label}Combined_R2.fastq.gz"

#done

#alignment
for dir in ficifolium/V1/F2toV1/*/
do
	Label=$(basename $dir)
    read1=$(find $Label -name "*R1*")
    read2=$(find $Label -name "*R2*")
    echo $read1
    echo $read2
    FileNameSam="${Label}.sam"
    FileNameBam="${Label}.bam"
    FileNameBamSorted="${Label}Sorted.bam"
    FileNameA="${Label}.fasta"
    FileNameQ="${Label}.fastq"
    Coverage="${Label}-Coverage.svg"
    BaiFile="${FileNameBam}.bai"
    VCFFile="${Label}.g.vcf.gz"
    VarBam="${Label}.variant.bam"

    #This sorts the cat'd zipped fastq files so BWA can properly align everything
    bash tools/bbmap/repair.sh in1=$read1 in2=$read2 out1="${Label}_Sorted_Combined_R1.fastq.gz" out2="${Label}_Sorted_Combined_R2.fastq.gz"
    mv "${Label}_Sorted_Combined_R1.fastq.gz" ${read1}
    mv "${Label}_Sorted_Combined_R2.fastq.gz" ${read2}

    #Alignment and sorting
    echo "BWA Using ${read1} and ${read2}"
    bwa mem -t $threads $RefLoc $read1 $read2 > $FileNameSam	#align reads to ref
    samtools view -@ $threads -b $FileNameSam > $FileNameBam	#convert sam to bam
    samtools sort -@ $threads $FileNameBam > $FileNameBamSorted	#sort bam
    mv $FileNameBamSorted $FileNameBam
    samtools index -b -@ $threads $FileNameBam

    #Coverage plot
    java -jar tools/jvarkit.jar wgscoverageplotter --include-contig-regex "Cf.*" --dimension 2500x800 -C 30 --clip -R $RefLoc $FileNameBam --percentile average > $Coverage

    #cleanup and move files to correct directory
    rm $FileNameSam
    mv $BaiFile "${dir}${BaiFile}"
    mv $FileNameBam "${dir}${FileNameBam}"
    mv $Coverage "${dir}${Coverage}"
done


