F2s were placed in individual folders containing forward and reverse reads. AlignAllReads.sh was called which recursively aligns each set of reads to the reference genome,
produces a BAM file, index, and coverage graph. ReheaderAndPileup.sh adjusts header information in each BAM file to reflect the folder the file is in
such that Freebayes can produce a VCF and properly associate BAM file data with the respective sample, and calls Freebayes-Parallel to produce the VCF. Finally
VCF2GWAS is passed this VCF file, and the phenotypic data to run the analyses. 
