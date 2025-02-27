#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --mem=92gb
#SBATCH --account=PAS1117
#SBATCH --job-name=dramv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1-47%10

assembly=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)

module use /fs/project/PAS1117/modulefiles
module load DRAM

# Annotate
fasta_input="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2-pass2/${assembly}/for-dramv/final-viral-combined-for-dramv.fa"
affi_input="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2-pass2/${assembly}/for-dramv/viral-affi-contigs-for-dramv.tab"

# Variables to pass to DRAMv annotate
opts="--skip_trnascan --threads 20 --min_contig_size 1000"
outDir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/DRAMv/${assembly}/DRAMv-annotate"
distillDir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/DRAMv/${assembly}/DRAMv-distill"

mkdir -p /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/DRAMv/${assembly}

time DRAM-v.py annotate -i $fasta_input -v $affi_input -o $outDir $opts

# Then summarize
time DRAM-v.py distill -i $outDir/annotations.tsv -o ${distillDir}
