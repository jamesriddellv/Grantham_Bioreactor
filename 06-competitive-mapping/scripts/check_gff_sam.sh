#!/bin/bash

cd ../01-combined-database

# Paths to your input file and output files
input_file="gff_headers.txt"
sam_file="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/STM_0716_E_M_E002_80FILTERED_NAMESORTED.sam"
included_file="included.txt"
excluded_file="excluded.txt"

# Clear the contents of the included and excluded files before starting
> $included_file
> $excluded_file

# Read each line from the input file
while read -r line; do
    # Check if the line exists in the SAM file
    if grep -qF "$line" "$sam_file"; then
        # If found, add the line to included.txt
        echo "$line" >> $included_file
    else
        # If not found, add the line to excluded.txt
        echo "$line" >> $excluded_file
    fi
done < "$input_file"
