#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=prodigal
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

module use /fs/project/PAS1117/modulefiles
module load prodigal/2.6.3

prodigal \
-i ../01-combined-database/filtered_clusters.fna \
-f gff \
-o ../01-combined-database/filtered_clusters.gff \
-a ../01-combined-database/filtered_clusters.faa \
-p meta

