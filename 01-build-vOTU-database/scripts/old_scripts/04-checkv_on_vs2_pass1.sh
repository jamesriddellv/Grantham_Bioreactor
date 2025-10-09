#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=checkv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out

######################
### load variables ###
######################

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

inFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2_pass1/${sample}/final-viral-combined.fa
outDir=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_vs2/${sample}
mkdir -p $outDir

checkvLoc="/users/PAS1117/osu9664/eMicro-Apps/CheckV-0.8.1.sif"

module load singularity
time $checkvLoc end_to_end $inFile $outDir -t 30

cat $outDir/proviruses.fna $outDir/viruses.fna > $outDir/combined.fna
