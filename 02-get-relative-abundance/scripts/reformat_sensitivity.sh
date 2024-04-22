#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=reformat_senstivitiy
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=1
#SBATCH --output=slurm-%x_%A_%a.out

module use /fs/project/PAS1117/modulefiles
module load singularityImages
module load samtools


sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaT.conf)

WORKING_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results
bt2_index_loc=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/bt2-index/checkv_10kb_or_high_quality_vOTUs
bam_dir=/fs/scratch/Sullivan_Lab/JamesR/checkv_10kb_or_high_quality_vOTUs/bam-files
bam_loc=${bam_dir}/${sample}.bam

echo  "reformatting with BBTools reformat.sh"
for i in $(seq 5 5 100)
do
    echo "filtering at ${i} ANI..."
    filtered_bam_loc=${bam_dir}/${sample}_${i}FILTERED.bam
    BBTools-38.97.sif reformat.sh -Xmx100g minidfilter=0.${i} in=${bam_loc} out=${filtered_bam_loc} pairedonly=t primaryonly=t overwrite=true
done
