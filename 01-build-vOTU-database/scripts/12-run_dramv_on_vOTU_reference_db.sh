#!/bin/bash
#SBATCH --time=01-00:00:00
#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --mem=92gb
#SBATCH --account=PAS1117
#SBATCH --job-name=dramv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

module use /fs/project/PAS1117/modulefiles
module load DRAM

# Annotate
fasta_input=../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/for-dramv/final-viral-combined-for-dramv.fa
affi_input=../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/for-dramv/viral-affi-contigs-for-dramv.tab

# Variables to pass to DRAMv annotate
opts="--skip_trnascan --threads 20 --min_contig_size 1000"
annoDir="../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/DRAMv/"
distillDir="../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/DRAMv/DRAMv-distill"

time DRAM-v.py annotate -i $fasta_input -v $affi_input -o $annoDir $opts

# Then summarize
# time DRAM-v.py distill -i $annoDir/annotations.tsv -o $distillDir
