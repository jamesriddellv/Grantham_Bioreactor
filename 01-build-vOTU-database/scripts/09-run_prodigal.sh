#!/bin/bash
#SBATCH --time=08:00:00
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
-i ../results/vOTU_clusters/checkv_5kb_or_high_quality_vOTUs.fna \
-f gff \
-o ../results/vOTU_clusters/checkv_5kb_or_high_quality_vOTUs.gff \
-a ../results/vOTU_clusters/checkv_5kb_or_high_quality_vOTUs.faa \
-p meta
