#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=bt2-build
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

module use /fs/project/PAS1117/modulefiles
module load bowtie2/2.4.1

bowtie2-build \
../../01-build-vOTU-database/results/vOTU_clusters/checkv_10kb_or_high_quality_vOTUs.fna \
../data/bt2-index/checkv_10kb_or_high_quality_vOTUs
