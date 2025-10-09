#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --account=PAS1117
#SBATCH --job-name=mvip
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


export CUDA_VISIBLE_DEVICES=-1
wDir=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/mvip-pipeline
metadata=../data/metadata_for_mvip.txt

### make metadata file for MVP pipeline ###

awk 'BEGIN {OFS="\t"; print "Sample_number", "Sample", "Assembly_Path", "Read_Path"} {split($0, a, "/"); print NR, a[length(a)], "../data/grantham_assemblies_5kb/" $0, "../data/grantham_assemblies_5kb/" $0}' ../assemblies_5kb.conf > ${metadata}
mkdir -p ${wDir}

### Run MVP pipeline ###

# mamba activate mvip

run MVP setup

mvip MVP_00_set_up_MVP \
-i ${wDir} \
-m ${metadata} \
--skip_check_errors \
--genomad_db_path /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/mvip-pipeline/00_DATABASES/genomad_db \
--checkv_db_path /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/mvip-pipeline/00_DATABASES/checkv-db-v1.5 \

# run MVP module 1

mvip MVP_01_run_genomad_checkv \
-i ${wDir} \
-m ${metadata} \
--min_seq_size 5000 \
--threads 30 \

# run MVP module 2

mvip MVP_02_filter_genomad_checkv \
-i ${wDir} \
-m ${metadata} \
--viral_min_genes 0 \
--host_viral_genes_ratio 1

# run MVP module 3

mvip MVP_03_do_clustering \
-i ${wDir} \
-m ${metadata} \
--threads 30

# run MVP module 6

mvip MVP_06_do_functional_annotation \
-i ${wDir} \
-m ${metadata} \
--DRAM \
--delete_files \
--threads 30
