#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH -n 16
#SBATCH --account=PAS1117
#SBATCH --job-name=bowtie2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=2-38%10
#SBATCH --output=slurm-%x_%A_%a.out

module use /fs/project/PAS1117/modulefiles
module load bowtie2/2.4.1
module load samtools/1.9



sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaT.conf)

R1_loc=trimmed_100M/${sample}_JGI_MT_R1_bbdtrimmed_100M.fastq
R2_loc=trimmed_100M/${sample}_JGI_MT_R2_bbdtrimmed_100M.fastq

WORKING_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results
SCRATCH_DIR=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor

echo "moving to the scratch directory..."

NUM_THREADS=16
bt2_index_loc=${WORKING_DIR}/bt2-index/vOTU_clusters
sam_loc=$TMPDIR/${sample}.sam

bam_dir=$WORKING_DIR/bam-files/${sample}

mkdir -p ${bam_dir}

bam_loc=${bam_dir}/${sample}.bam

cd $TMPDIR

echo "running bowtie2..."
bowtie2 --threads $NUM_THREADS -x $bt2_index_loc -1 ${SCRATCH_DIR}/metaT/$R1_loc -2 ${SCRATCH_DIR}/metaT/$R2_loc -S $sam_loc
samtools sort ${sam_loc} -o ${bam_loc}
