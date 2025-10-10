#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH -n 16
#SBATCH --account=PAS1117
#SBATCH --job-name=bt2-build
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


inFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/prefix_MAGs.fasta
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/data/bt2-index

cat /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/MAGs/prefix/* > $inFile

module load bowtie2

bowtie2-build \
${inFile} \
${outDir}/prefix_MAGs \
--threads 16 \
--seed 42 \
--large-index

rm $inFile
