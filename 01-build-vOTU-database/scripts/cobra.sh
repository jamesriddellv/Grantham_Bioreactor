#!/bin/bash
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -A PAS1117
#SBATCH -J cobra


module use /fs/ess/PAS1117/modulefiles
module load COBRA/1.1.54

Workdir=/fs/ess/PAS1117/Grantham/20230407_Grantham_incubations/data_processing/COBRA
cd $Workdir

cobra.py \
-f ../results/vOTU_clusters/vOTU_clusters_5kb_cutoff.fna \
-q ../results/vOTU_clusters/vOTU_clusters_5kb_cutoff.fna \
-c ./01_read_mapping/STM_0716_E_M_E061/STM_0716_E_M_E061_meancov.txt \
-m ./01_read_mapping/STM_0716_E_M_E061/bam_files/STM_0716_E_M_E061_megahit_final.contigs_1000.fa.STM_0716_E_M_E061_R1_trimmed.fastq.gz.bam \
-a megahit \
-mink 31 \
-maxk 121 \
-o STM_0716_E_M_E061.COBRA.out
