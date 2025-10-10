#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -n 40
#SBATCH --account=PAS1117
#SBATCH --job-name=bowtie2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=1-10
#SBATCH --output=slurm-%x_%A_%a.out

module use /fs/project/PAS1117/modulefiles
module load singularityImages
module load bowtie2
module load samtools

bowtie2 --version
samtools --version

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaG.conf)

R1_loc=${sample}_R1_trimmed.fastq
R2_loc=${sample}_R2_trimmed.fastq
NUM_THREADS=40
bt2_index_loc=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/bt2-index/prefix_MAGs
sam_loc=${sample}.sam
bam_dir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaG/bam
bam_loc=${bam_dir}/${sample}.bam
filtered_bam_loc=${bam_dir}/${sample}_97FILTERED.bam
filtered_sorted_bam_loc=${bam_dir}/${sample}_97FILTERED_NAMESORTED

mkdir -p ${bam_dir}

cd /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/050724_GranthamMetaG

zcat ${R1_loc}.gz > ${R1_loc}
zcat ${R2_loc}.gz > ${R2_loc}

echo ${R1_loc}
echo ${R2_loc}

echo "running bowtie2..."
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p $NUM_THREADS -x $bt2_index_loc -1 $R1_loc -2 $R2_loc -S $sam_loc
samtools view -@ ${NUM_THREADS} -bS ${sam_loc} > ${bam_loc}

echo  "reformatting with BBTools reformat.sh"
echo "filtering at 97% ANI..."

BBTools-38.97.sif reformat.sh -Xmx100g minidfilter=0.97 in=${bam_loc} out=${filtered_bam_loc} pairedonly=t primaryonly=t overwrite=true
samtools sort ${filtered_bam_loc} -o ${filtered_sorted_bam_loc}.bam -@ ${NUM_THREADS}

echo "read mapping complete."

rm ${R1_loc} ${R2_loc}
