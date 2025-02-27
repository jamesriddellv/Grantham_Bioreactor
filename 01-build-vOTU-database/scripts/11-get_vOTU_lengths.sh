#!/bin/bash

# Hardcoded input and output file paths
input_fasta="../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs/genomad_vs2sop_sun2024_emerson2018_combined.fa.self-blastn.clusters.fna"
output_tsv="../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs/vOTU_contig_lengths.tsv"

# Check if input file exists
if [ ! -f "$input_fasta" ]; then
    echo "Error: Input FASTA file '$input_fasta' not found."
    exit 1
fi

# Process the FASTA file and generate the TSV file
awk '
    BEGIN {
        FS = "\n";
        RS = ">";
        OFS = "\t";
    }
    NR > 1 {
        header = $1;
        seq = "";
        for (i = 2; i <= NF; i++) {
            seq = seq $i;
        }
        print header, length(seq);
    }
' "$input_fasta" > "$output_tsv"

echo "TSV file with sequence lengths generated: $output_tsv"

