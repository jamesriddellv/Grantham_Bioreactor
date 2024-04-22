#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=htseq-count
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


# run htseq-count
htseq-count -a 0 \
-t CDS \
-i ID \
-n 28 \
--stranded=no \
-c ../03-mapped-reads-per-gene/htseq_vOTUs_100M_80FILTERED_UNSTRANDED.tsv \
../02-mapping-results/bam-files/*80FILTERED_NAMESORTED.bam \
../01-combined-database/filtered_clusters_host_genes.gff

bash remove_zero_rows_htseq_table.sh ../03-mapped-reads-per-gene/htseq_vOTUs_100M_80FILTERED_UNSTRANDED.tsv ../03-mapped-reads-per-gene/htseq_vOTUs_100M_80FILTERED_UNSTRANDED_no0s.tsv
