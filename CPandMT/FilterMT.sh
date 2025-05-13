#!/bin/bash

#takes original freebayes vcf and filters it down for further processing
#required settings 
filename="RawMT"
numsamples=6

#tunable parameters
MinDepthPerSampleAverage=100 # on average, at least how many reads should be present per individual
MaxDepthPerSampleAverage=3000 # on average, what is the upper limit for reads per individual
MinQual=998 #minimum site quality
minGQ=30 # what is the minimum quality for a call to be kept
MinDP=50 # minimum depth per sample for that call to be retained
MaxMissing=0.5 # minimum percent of individuals called to keep the locus

#math section
mindp=$(( numsamples * MinDepthPerSampleAverage ))
maxdp=$(( numsamples * MaxDepthPerSampleAverage ))

#Run the first filtering step with vcftools
bcftools stats $filename.vcf > $filename.vcf.stats
vcftools --vcf $filename.vcf --minDP $MinDP --recode --recode-INFO-all --out $filename.GQ.vcf
mv $filename.GQ.vcf.recode.vcf $filename.GQ.vcf #move the output file because vcftools renames it
bcftools stats $filename.GQ.vcf > $filename.GQ.vcf.stats

#Filter on depths and quality
bcftools filter -o $filename.DPGQ.vcf -e 'INFO/DP < '$mindp' || INFO/DP > '$maxdp' || QUAL < '$MinQual $filename.GQ.vcf
bcftools stats $filename.DPGQ.vcf > $filename.DPGQ.vcf.stats

#remove any sites with allele balance below 25%, above 75% or about 0. Following lines pulled from http://www.ddocent.com/filtering/
vcffilter -s -f "AB > 0.25 & AB < 0.75 | AB < 0.01" $filename.DPGQ.vcf > $filename.DPGQAB.vcf
bcftools stats $filename.DPGQAB.vcf > $filename.DPGQAB.vcf.stats

#filter sites from both strands
#vcffilter -f "SAF / SAR > 100 & SRF / SRR > 100 | SAR / SAF > 100 & SRR / SRF > 100" -s $filename.DPGQAB.vcf > $filename.DPGQAB.strand.vcf # failed
#bcftools stats $filename.DPGQAB.strand.vcf > $filename.DPGQAB.strand.vcf.stats

#filter on quality of reads supporing both the ref and alt alleles
vcffilter -f "PAIRED > 0.05 & PAIREDR > 0.05 & PAIREDR / PAIRED < 1.75 & PAIREDR / PAIRED > 0.25 | PAIRED < 0.05 & PAIREDR < 0.05" -s $filename.DPGQAB.vcf > $filename.DPGQAB.paired.vcf

#use the following line for RADseq and GBS - filter if there is a large difference between sampled reference and alt alleles
#vcffilter -f "MQM / MQMR > 0.9 & MQM / MQMR < 1.05" in.vcf > out.vcf

#print stats of output file
bcftools stats $filename.DPGQAB.paired.vcf > $filename.DPGQAB.paired.vcf.stats

#remove unused files
rm $filename.GQ.vcf
rm $filename.DPGQ.vcf
rm $filename.DPGQAB.vcf

