#!/bin/bash

# Input and output file names
input_file="60716-CDS-prot.fasta"
output_file="60716-CDS-prot-reformat.fasta"

# Process the input file
awk '
    BEGIN { FS="\\|\\|" } 
    /^>/ { 
        split($0, a, "\\|\\|")
        print ">" a[5] 
    } 
    !/^>/ { print }' $input_file > $output_file

echo "Reformatting complete. Output written to $output_file"
