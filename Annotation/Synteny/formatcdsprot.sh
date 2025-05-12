#!/bin/bash

# Input and output file paths
input_file="protein.faa"
output_file="10.1.protein.cleaned.faa"


awk '/^>/{sub(/ .*/, ""); print; next} {print}' "$input_file" > "$output_file"