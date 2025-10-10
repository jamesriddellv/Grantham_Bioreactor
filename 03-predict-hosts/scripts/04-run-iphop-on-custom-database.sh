#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --account=PAS1117
#SBATCH --job-name=iphop
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module use /fs/project/PAS1117/modulefiles
module load iPHoP/1.3.3

inFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/03-predict-hosts/results/iphop_custom
dbDir=../data/prefix/iphop_Aug_2023_pub_prefix_grantham_rw
THREADS=30

mkdir -p ${outDir}

iphop predict \
--fa_file ${inFile} \
--db_dir ${dbDir} \
--out_dir ${outDir} \
-t ${THREADS}

