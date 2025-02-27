#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=cluster_viruses
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=%x_%j.out

### combine vs2_sop, genomad, stordalen, emerson_2018
genomad=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/all_viral_contigs.fna
vs2sop=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2-sop-screen/vs2sop-final-viral-scored.fa
sun2024=../data/stordalen_vOTUs_symlink.fa
emerson2018=../data/emerson_2018_vOTUs.fa

inFile=../results/viral_contigs/genomad_vs2sop_sun2024_emerson2018_combined.fa
outDir=../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs

cat ${genomad} ${vs2sop} ${sun2024} ${emerson2018} > ${combinedFile}

mkdir -p $outDir

### run checkv on combined viral contigs
module use /fs/project/PAS1117/modulefiles
module load singularityImages

CheckV-0.8.1-ClusterONLY.sif \
-i $inFile \
-t 30 \
-o $outDir \
--min-ani 95 \
--min-qcov 0 \
--min-tcov 85
