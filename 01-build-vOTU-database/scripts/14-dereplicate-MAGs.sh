#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=galah
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

cd ../data

galah cluster \
-f MAGs/* \
--output-representative-fasta-directory MAGs_galah \
--precluster-ani 90 \
--ani 97 \
--precluster-method finch \
--threads 28
