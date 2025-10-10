#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH -A PAS1117
#SBATCH -J coverm_metaG
#SBATCH --output=slurm-%x_%j.out
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL

mode=trimmed_mean
ANI=90
ALIGNED_PERC=75
MIN_COVERED_FRACTION=10

BAM_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaG/bam
OUT_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaG/coverm_trimmed_mean

mkdir -p ${OUT_DIR}

module load singularity/current

while read -r sample; do
    singularity run /users/PAS1117/osu9664/eMicro-Apps/CoverM-0.6.1.sif \
    contig \
    --bam-files ${BAM_DIR}/${sample}_90FILTERED_SORTED.bam \
    --output-format dense \
    --min-read-percent-identity ${ANI} \
    --min-read-aligned-percent ${ALIGNED_PERC} \
    --min-covered-fraction ${MIN_COVERED_FRACTION} \
    -m ${mode} \
    -t 30 > ${OUT_DIR}/${sample}_${mode}_${ANI}_${ALIGNED_PERC}_${MIN_COVERED_FRACTION}.txt
done < ./metaG.conf
