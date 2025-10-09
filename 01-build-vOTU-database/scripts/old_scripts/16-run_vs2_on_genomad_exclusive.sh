#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH -n 28
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass1
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
inFile="../results/vOTU_clusters/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vs2_exclusive_vOTUs.fna"
outDir="../results/vOTU_clusters/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/genomad_exlusive_vs2_classification"
minLength=5000
opts="--keep-original-seq --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-score 0.01 -j 28"

mkdir -p $outDir

module load singularity
time $vs2Loc run \
-i ${inFile} \
-w $outDir \
--min-length $minLength \
$opts all --rerun-incomplete
