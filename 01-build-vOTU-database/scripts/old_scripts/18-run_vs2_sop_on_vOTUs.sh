#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
inFile="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTU_clusters/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vs2_sop_genomad_stordalen_emerson2018_combined.fa.self-blastn.clusters.fna"
outDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vs2_no_cutoff"
minLength=5000
opts="--keep-original-seq --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-score 0.0 -j 28"

module load singularity
time $vs2Loc run \
-i $inFile \
-w $outDir \
--min-length $minLength \
$opts all

wait
