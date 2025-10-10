#!/bin/bash
#SBATCH -t 00:10:00
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -A PAS1117
#SBATCH -J coverm
#SBATCH --array=1-10
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaG.conf)

mode=mean
ANI=97
ALIGNED_PERC=75
MIN_COVERED_FRACTION=75

BAM_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaG/bam
OUT_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaG/coverm

mkdir -p ${OUT_DIR}

module load singularity/current


singularity run /users/PAS1117/osu9664/eMicro-Apps/CoverM-0.6.1.sif \
genome \
--bam-files ${BAM_DIR}/${sample}_97FILTERED_SORTED.bam \
--genome-fasta-directory ../../01-build-vOTU-database/data/MAGs/prefix \
-x "fa" \
--output-format dense \
--min-read-percent-identity ${ANI} \
--min-read-aligned-percent ${ALIGNED_PERC} \
--min-covered-fraction ${MIN_COVERED_FRACTION} \
-m ${mode} \
-t 40 > ${OUT_DIR}/${sample}_${mode}_${ANI}_${ALIGNED_PERC}_${MIN_COVERED_FRACTION}.txt
