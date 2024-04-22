#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=htseq-count
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


outFile="../results/htseq_5kb_vOTUs_100M_85FILTERED_REVSTRANDED.tsv"
annoFile="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTU_clusters/checkv_5kb_or_high_quality_vOTUs.gff"
bam_dir="/fs/scratch/Sullivan_Lab/JamesR/checkv_5kb_or_high_quality_vOTUs/bam-files"

no_zeroes="../results/htseq_5kb_vOTUs_100M_85FILTERED_REVSTRANDED_no0s.tsv"


# run htseq-count
htseq-count -a 0 \
-t CDS \
-i ID \
-n 28 \
--stranded=reverse \
-c $outFile \
${bam_dir}/*85FILTERED_NAMESORTED.bam \
$annoFile

bash remove_zero_rows_htseq_table.sh $outFile $no_zeroes
