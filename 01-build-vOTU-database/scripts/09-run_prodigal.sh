#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=prodigal
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out
module use /fs/project/PAS1117/modulefiles
module load prodigal/2.6.3
prodigal \
-i ../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vs2_sop_genomad_stordalen_emerson2018_combined.fa.self-blastn.clusters.fna \
-f gff \
-o ../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/final-viral-combined.gff \
-a ../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/final-viral-combined.faa \
-d ../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/final-viral-combined-genes.fna \
-p meta
