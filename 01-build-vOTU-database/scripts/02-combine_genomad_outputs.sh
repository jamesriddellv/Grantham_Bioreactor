# combines genomad virus predictions from each assembly

# fasta
cat /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/*/*_summary/*_virus.fna > /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/all_viral_contigs.fna

# summary stats
cat /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/*/*_summary/*_virus_summary.tsv > /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/all_viral_contigs_summary.tsv
