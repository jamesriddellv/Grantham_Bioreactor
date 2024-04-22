cd ../01-combined-database

# Filter BLAST results based on criteria (>95% ANI and >85% coverage)
awk '$3 >= 95 && $11 >= 85' vOTU_against_host_db_results.txt | cut -f1 | sort -u > filtered_virus_ids.txt

# Use SeqKit to filter virus sequences
seqkit grep -v -f filtered_virus_ids.txt all_seqs.fna.self-blastn.clusters.fna > filtered_clusters.fna
