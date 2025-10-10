#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --account=PAS1117
#SBATCH --job-name=bt2-build
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


inFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/bt2-index
mkdir -p ${outDir}

module load bowtie2

bowtie2-build --version

bowtie2-build \
${inFile} \
${outDir}/combined_manual_filtered \
--threads 30 \
--seed 42
