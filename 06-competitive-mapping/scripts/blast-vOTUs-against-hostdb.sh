#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH -n 40
#SBATCH --account=PAS1117
#SBATCH --job-name=blastvOTU
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

cd ../01-combined-database

blastn \
-query all_seqs.fna.self-blastn.clusters.fna \
-db hosts_db \
-outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovs" \
-out vOTU_against_host_db_results.txt \
-num_threads 40
