#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH -n 10
#SBATCH --account=PAS1117
#SBATCH --job-name=VS2-SOP-screen
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=%x_%A_%a.out
#SBATCH --array=1-47

assembly=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

vs2_genomes="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_vs2/${sample}/combined.fna"
vs2_final_score="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2_pass1/${sample}/final-viral-score.tsv"
amg_summary="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/DRAMv/${assembly}/DRAMv-distill/amg_summary.tsv"
checkv_contamination="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_vs2/${sample}/contamination.tsv"

output_dir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2-sop-screen/${assembly}"

mkdir -p $output_dir

python /users/PAS1117/osu9664/eMicro-Apps/Process-VS2_and_DRAMv.py \
--vs2-scores $vs2_final_score \
--checkv-contam $checkv_contamination \
--dramv-amg $amg_summary \
--vs2-genomes $vs2_genomes \
--output-dir $output_dir \
--drop-manual
