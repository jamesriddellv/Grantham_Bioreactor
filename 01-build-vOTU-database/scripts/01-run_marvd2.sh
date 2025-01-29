#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

inFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/April2023_assemblies.fa
# inFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/October2023_assemblies.fa
# inFile=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/MAGs_assemblies.fa

outDir=../results/viral_id_outputs/marvd2/April2023
mkdir -p $outDir

module use /fs/project/PAS1117/modulefiles
module load MArVD2/0.11.8
