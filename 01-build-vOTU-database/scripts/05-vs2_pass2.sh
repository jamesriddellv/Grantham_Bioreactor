#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1-47%10

assembly=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

cd /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
inFile="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_vs2/${sample}/combined.fna"
outDir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2-pass2/${assembly}"
opts="--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -j 30 --include-groups dsDNAphage,ssDNA --min-score 0.01"

module load singularity
time $vs2Loc run \
-i ${inFile} \
-w $outDir \
$opts all

wait
