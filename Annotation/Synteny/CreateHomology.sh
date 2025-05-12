#!/bin/bash

input_blast="cfxcq.blast"
output_homology="cfcq.homology"

# Extract gene pairs and bit scores
awk '{print $1, $2, $12}' $input_blast > $output_homology

echo "Homology file created: $output_homology"