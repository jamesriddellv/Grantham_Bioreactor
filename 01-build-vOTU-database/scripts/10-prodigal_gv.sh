#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=prodigal
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


# mamba activate prodigal-gv

genomeFile=../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs/genomad_vs2sop_sun2024_emerson2018_combined.fa.self-blastn.clusters.fna
outDir=../results/genomad_vs2sop_sun2024_emerson2018_combined_vOTUs

python3 parallel-prodigal-gv.py \
-t 30 \
-q \
-p meta \
-i ${genomeFile} \
-f gff \
-o ${outDir}/vOTU_genes.gff \
-a ${outDir}/vOTU_genes.faa \
-d ${outDir}/vOTU_genes.fna \
