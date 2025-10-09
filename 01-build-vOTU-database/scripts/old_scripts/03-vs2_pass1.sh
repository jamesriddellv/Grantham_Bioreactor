#!/bin/bash
#SBATCH --time=04-00:00:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass1
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1

######################
### load variables ###
######################

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

inFile=../data/grantham_assemblies_5kb/${sample}
outDir=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2_pass1/${sample}
mkdir -p $outDir

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"

minLength=5000
opts="--keep-original-seq --include-groups dsDNAphage,ssDNA --min-score 0.5 -j 30"

#############################
### run VirSorter2 Pass 1 ###
#############################
module load singularity
time $vs2Loc run \
-i ${inFile} \
-w ${outDir} \
--min-length ${minLength} \
$opts all --rerun-incomplete
