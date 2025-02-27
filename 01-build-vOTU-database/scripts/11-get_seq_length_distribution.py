from Bio import SeqIO
import numpy as np

def get_sequence_lengths(fasta_file):
    lengths = [len(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]
    return lengths

def print_summary_statistics(lengths):
    print(f"Total sequences: {len(lengths)}")
    print(f"Minimum length: {np.min(lengths)}")
    print(f"Maximum length: {np.max(lengths)}")
    print(f"25th percentile: {np.percentile(lengths, 25)}")
    print(f"50th percentile (median): {np.percentile(lengths, 50)}")
    print(f"75th percentile: {np.percentile(lengths, 75)}")

if __name__ == "__main__":
   fasta_file_path = "../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs/genomad_vs2sop_sun2024_emerson2018_combined.fa.self-blastn.clusters.fna"
   sequence_lengths = get_sequence_lengths(fasta_file_path)
   print_summary_statistics(sequence_lengths)
