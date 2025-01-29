#!/bin/bash

# Directory containing the fasta files
input_dir="../data/MAGs"

mkdir -p ../data/MAGs/prefix

# Loop through each .fa file in the directory
for file in "$input_dir"/*.fa; do
    # Get the base name of the file
    base_name=$(basename "$file")
    
    # Create the new filename with "Add_" prefix
    new_file="${input_dir}/prefix/Add_${base_name}"
    
    # Add "Add_" to each header and save to the new file
    awk '/^>/ {print "Add_" $0; next} {print}' "$file" > "$new_file"
done
