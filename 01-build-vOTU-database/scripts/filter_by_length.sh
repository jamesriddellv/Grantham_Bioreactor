# purpose is to filter all sequences above 3kb and then extend these with COBRA.
seqkit seq --min-len 3000 -w 0 ../results/vOTU_clusters/all_seqs.fna.self-blastn.clusters.fna > ../results/vOTU_clusters/vOTU_clusters_5kb_cutoff.fna
