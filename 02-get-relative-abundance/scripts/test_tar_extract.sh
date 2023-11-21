#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --account=PAS1117
#SBATCH --job-name=bowtie2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=1-2
#SBATCH --output=slurm-%x_%A_%a.out


sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaT.conf)

R1_loc=trimmed_100M/${sample}_JGI_MT_R1_bbdtrimmed_100M.fastq
R2_loc=trimmed_100M/${sample}_JGI_MT_R2_bbdtrimmed_100M.fastq

WORKING_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/
SCRATCH_DIR=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/

echo "moving to the scratch directory..."

cd ${SCRATCH_DIR}/metaT

# Extract R1 and R2
echo "extracting ${sample} R1 and R2..."

tar -xvf ${SCRATCH_DIR}/metaT/trimmed_100M.tar.gz ${R1_loc}
