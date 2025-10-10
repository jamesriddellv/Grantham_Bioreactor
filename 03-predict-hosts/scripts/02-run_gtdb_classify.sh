#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=gtdb_classify
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

dataDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/MAGs/prefix"
outDir="../data/gtdb_classify_wf"
threads=40

gtdbtk classify_wf \
--genome_dir ${dataDir}/ \
--out_dir ${outDir} \
--cpus ${threads} \
--pplacer_cpus ${threads} \
--extension fasta \
--skip_ani_screen
