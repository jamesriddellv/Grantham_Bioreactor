#!/bin/bash

# bash 10-get_vOTU_lengths.sh ../results/vOTU_clusters/vs2_sop_viral_contigs_stordalen_vOTUs_combined/vs2_sop_viral_contigs_stordalen_vOTUs_combined.fa.self-blastn.clusters.fna ../results/vOTU_clusters/vs2_sop_viral_contigs_stordalen_vOTUs_combined/vs2_sop_stordalen_vOTUs_contig_lengths.tsv

# Input FASTA file
input_fasta="$1"

# Output TSV file
output_tsv="$2"

# Check if input and output files are provided
if [ -z "$input_fasta" ] || [ -z "$output_tsv" ]; then
    echo "Usage: $0 input_fasta output_tsv"
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

