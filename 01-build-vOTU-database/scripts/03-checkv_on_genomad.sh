#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=checkv-on-genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1-47%10

######################
### load variables ###
######################

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

sample_base=$(basename "${sample}" .fa)

inFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/${sample}/${sample_base}_summary/${sample_base}_virus.fna
outDir=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad/${sample}
mkdir -p $outDir

checkv -h

time checkv end_to_end $inFile $outDir -t 30

# Handle situations where CheckV identifies additional proviruses

echo "number of additional proviruses detected that need to be re-classified by genomad"
grep ">" $outDir/proviruses.fna | wc -l

# cat $outDir/proviruses.fna $outDir/viruses.fna > $outDir/combined.fna
