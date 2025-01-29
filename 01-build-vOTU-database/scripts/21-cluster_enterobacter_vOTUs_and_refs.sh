#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=cluster_viruses
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=%x_%j.out

# default

inFile=../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/enterobacter_vOTU_and_vcontact3_cluster_references.fna
outDir=../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/enterobacter_vOTU_and_vcontact3_cluster_references

mkdir -p $outDir

module use /fs/project/PAS1117/modulefiles
module load singularityImages

CheckV-0.8.1-ClusterONLY.sif \
-i $inFile \
-t 48 \
-o $outDir \
--min-ani 95 \
--min-qcov 0 \
--min-tcov 85
