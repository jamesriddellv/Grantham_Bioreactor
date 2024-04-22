#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=bowtie2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --array=1
#SBATCH --output=slurm-%x_%A_%a.out

# This script is a sensitivity analysis of reformat.sh. Bridget McGivern had used 97% ANI cutoffs for the reads mapping back,
# but that was for microbes, and here we are doing viruses. The goal is to be able to compare viral transcript abundance to
# microbial transcript abundance, but since our viruses are clustered at 95% ANI 85% qcov, they may have valid reads mapping
# at 80% ANI because a gene not represented fully from the cluster representative is still in that viral cluster.

# Here we'll test 0-100% cutoffs with increments of 10 for one or two samples and see how the abundance changes, and see if
# there's any clear "fall off." I can go over it with Matt too.

module use /fs/project/PAS1117/modulefiles
module load bowtie2/2.4.1
module load samtools/1.16.1



sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ./metaT.conf)


R1_loc=trimmed_100M/${sample}_JGI_MT_R1_bbdtrimmed_100M.fastq
R2_loc=trimmed_100M/${sample}_JGI_MT_R2_bbdtrimmed_100M.fastq


WORKING_DIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/06-competitive-mapping/02-mapping-results


NUM_THREADS=28
bt2_index_loc=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/06-competitive-mapping/01-combined-database/bt2-index/vOTU-host-combined
sam_loc=$TMPDIR/${sample}.sam

bam_dir=$WORKING_DIR/sensitivity-analysis/${sample}

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


echo  "reformatting with BBTools reformat and filtering to 10-100% identity every 10..."
module use /fs/project/PAS1117/modulefiles
module load singularityImages

for i in $(seq 0 10 100); do
    filtered_bam_loc=${bam_dir}/${sample}_${i}FILTERED.bam
    BBTools-38.97.sif reformat.sh -Xmx100g idfilter=0.${i} in=${bam_loc} out=${filtered_bam_loc} pairedonly=t primaryonly=t overwrite=true
    echo "sorting 0.${i} filtered bam..."
    filtered_sorted_bam_loc=${bam_dir}/${sample}_${i}FILTERED_NAMESORTED
    samtools sort -n -T ${filtered_sorted_bam_loc} -o ${filtered_sorted_bam_loc}.bam ${filtered_bam_loc} -@ ${NUM_THREADS}
    rm ${filtered_bam_loc}
done
echo "read mapping complete."

