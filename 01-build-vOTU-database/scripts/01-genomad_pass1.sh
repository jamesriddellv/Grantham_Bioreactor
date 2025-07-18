#!/bin/bash

#########################################
### Specify slurm batch job variables ###
#########################################

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1-49%10

######################
### load variables ###
######################

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies_5kb.conf)

inFile=../data/grantham_assemblies_5kb/${sample}
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/${sample}
dbDir=../data/genomad_db

mkdir -p $outDir

###################
### run genomad ###
###################

genomad end-to-end --cleanup --enable-score-calibration --threads 30 $inFile $outDir $dbDir

