#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --account=PAS1117
#SBATCH --job-name=dramv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

module use /fs/project/PAS1117/modulefiles
module load DRAM
# Annotate
fasta_input="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/mvip/06_FUNCTIONAL_ANNOTATION/06_DRAM_V/MVP_03_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Sequences_DRAM_Input.fna"
affi_input="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/mvip/06_FUNCTIONAL_ANNOTATION/06_DRAM_V/MVP_06_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Gene_Annotation_GENOMAD_DRAM_Input.tsv"

# Variables to pass to DRAMv annotate
opts="--skip_trnascan --threads 30 --min_contig_size 1000"
outDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/DRAMv"
annoDir="${outDir}/DRAMv-annotate"
distillDir="${outDir}/DRAMv-distill"

mkdir -p ${outDir}

time DRAM-v.py annotate -i $fasta_input -v $affi_input -o $annoDir $opts

# Then summarize
time DRAM-v.py distill -i $annoDir/annotations.tsv -o $distillDir
