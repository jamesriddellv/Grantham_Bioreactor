#!/bin/bash

#########################################
### Specify slurm batch job variables ###
#########################################

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1-47%10

######################
### load variables ###
######################

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

inFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad/${sample}/proviruses.fna
outDir=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad_pass2/${sample}
dbDir=../data/genomad_db

mkdir -p $outDir

###################
### run genomad ###
###################

genomad end-to-end --cleanup --enable-score-calibration --threads 30 $inFile $outDir $dbDir

