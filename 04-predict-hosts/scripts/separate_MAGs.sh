#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -n 8
#SBATCH --account=PAS1117
#SBATCH --job-name=separate_mags
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL


scaffolds="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/scaffolds_1kb_renamed/scaffolds.fa"
MagList="../data/Grantham_MAGs.list"

cat ${MagList} | while read line; do seqkit grep -w 0 -r -p ^${line}_ ${scaffolds} -o /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/MAGs/${line}.fa; done
