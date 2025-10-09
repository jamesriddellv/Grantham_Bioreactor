#!/bin/bash

#########################################
### Specify slurm batch job variables ###
#########################################

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=genomad_pass2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

######################
### load variables ###
######################

inFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined.fna.self-blastn.clusters.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/genomad_annotate
dbDir=../data/genomad_db

mkdir -p $outDir

###################
### run genomad ###
###################

genomad annotate $inFile $outDir $dbDir
