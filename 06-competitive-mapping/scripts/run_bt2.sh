#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=bowtie2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=1-38%10
#SBATCH --output=slurm-%x_%A_%a.out

module use /fs/project/PAS1117/modulefiles
module load bowtie2/2.4.1
module load samtools/1.16.1



sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaT.conf)


R1_loc=trimmed_100M/${sample}_JGI_MT_R1_bbdtrimmed_100M.fastq
R2_loc=trimmed_100M/${sample}_JGI_MT_R2_bbdtrimmed_100M.fastq


WORKING_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/06-competitive-mapping/02-mapping-results


NUM_THREADS=28
bt2_index_loc=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/06-competitive-mapping/scripts/filtered_clusters
sam_loc=$TMPDIR/${sample}.sam

bam_dir=$WORKING_DIR/filtered_clusters-bam-files

mkdir -p ${bam_dir}

bam_loc=${bam_dir}/${sample}.bam

cd $TMPDIR

echo "extracting metaT reads from tar archive to temporary directory..."

tar -xvf \
/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/06-competitive-mapping/01-combined-database/trimmed_100M.tar.gz \
${R1_loc}.gz \
${R2_loc}.gz \
-C ./

gunzip ${R1_loc}
gunzip ${R2_loc}

echo "running bowtie2..."
bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p $NUM_THREADS -x $bt2_index_loc -1 $R1_loc -2 $R2_loc -S $sam_loc

samtools view -@ ${NUM_THREADS} -bS ${sam_loc} > ${bam_loc}

filtered_bam_loc=${bam_dir}/${sample}_80FILTERED.bam

echo  "reformatting with BBTools reformat and filtering to 80% identity..."
module use /fs/project/PAS1117/modulefiles
module load singularityImages
BBTools-38.97.sif reformat.sh -Xmx100g idfilter=0.80 in=${bam_loc} out=${filtered_bam_loc} pairedonly=t primaryonly=t overwrite=true

echo "sorting filtered bam..."
filtered_sorted_bam_loc=${bam_dir}/${sample}_80FILTERED_NAMESORTED

samtools sort -n -T ${filtered_sorted_bam_loc} -o ${filtered_sorted_bam_loc}.bam ${filtered_bam_loc} -@ ${NUM_THREADS}

echo "sorting unfiltered bam..."

# also sort the unfiltered
samtools sort -n ${bam_loc} -o ${bam_dir}/${sample}_NAMESORTED.bam -@ ${NUM_THREADS}

echo "read mapping complete."
