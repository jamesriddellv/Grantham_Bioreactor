#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -n 8
#SBATCH --account=PAS1117
#SBATCH --job-name=prodigal
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out



# mamba activate prodigal-gv

genomeFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/data/prefix_MAGs.fasta
outDir=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/MAGs/prodigal

cat /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/MAGs/prefix/* > $genomeFile

mkdir -p ${outDir}

module use /fs/project/PAS1117/modulefiles
module load prodigal/2.6.3

prodigal -q \
-i ${genomeFile} \
-f gff \
-o ${outDir}/prefix_MAGs.gff \
-a ${outDir}/prefix_MAGs.faa \
-d ${outDir}/prefix_MAGs.fna \
