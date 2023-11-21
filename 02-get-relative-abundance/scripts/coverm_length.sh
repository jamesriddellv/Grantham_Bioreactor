#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 28
#SBATCH -A PAS1117
#SBATCH -J coverm
#SBATCH --array=1
#SBATCH --output=slurm-outputs/slurm-%x_%A_%a.out
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaT.conf)

WORKING_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/

module load singularity/current
singularity run /users/PAS1117/osu9664/eMicro-Apps/CoverM-0.6.1.sif \
contig \
--bam-files ${WORKING_DIR}/results/bam-files/${sample}/${sample}.bam \
--output-format dense \
-m length \
-t 28 > ${WORKING_DIR}/results/bam-files/vOTU_length.txt
