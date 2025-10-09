#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --account=PAS1117
#SBATCH --job-name=prodigal
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out



# mamba activate prodigal-gv

genomeFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/prodigal-gv

mkdir -p ${outDir}

python3 parallel-prodigal-gv.py \
-t 30 \
-q \
-i ${genomeFile} \
-f gff \
-o ${outDir}/combined_manual_filtered.gff \
-a ${outDir}/combined_manual_filtered.faa \
-d ${outDir}/combined_manual_filtered.fna \

