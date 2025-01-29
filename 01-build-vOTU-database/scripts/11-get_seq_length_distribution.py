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
#   fasta_file_path = "/fs/ess/PAS1117/Grantham/20230407_Grantham_incubations/VS2-SOP-screen/final-viral-scored-all.fa"  # Replace with the path to your FASTA file
#   fasta_file_path = "../data/2023-10-31-assemblies-final-viral-scored.fa" 
#   fasta_file_path = "/fs/ess/PAS1117/apratama/isogenie/Isogenie_final_viral_contigs/FINAL_5051_votus.fasta"
   fasta_file_path = "../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/present_metaT_vOTUs_post_manual_curation_1x.fna"
   sequence_lengths = get_sequence_lengths(fasta_file_path)
   print_summary_statistics(sequence_lengths)
