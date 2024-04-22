#!/bin/bash

if [ $# -lt 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Assign input file and output file from arguments
input_file="$1"
output_file="$2"

awk '{
    # Initialize a flag to determine if the row contains only zeros in all but the first column
    all_zeros = 1
    # Start checking columns from the second column to the end
    for (i = 2; i <= NF; i++) {
        # If a non-zero value is found, set the flag to 0
        if ($i != 0) {
            all_zeros = 0
            break
        }
    }
    # Print the line only if all_zeros flag is not set
    if (!all_zeros) {
        print
    }
}' "$input_file" > "$output_file"
