#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=vs2_pass2
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
inFile="../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vOTUs_no_detection_cutoff_post_manual_curation_simplified_headers.fna"
outDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/prep-for-dramv-test/"
minLength=0
opts="--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -j 40 --include-groups dsDNAphage --min-score 0.0"

module load singularity
time $vs2Loc run \
-i $inFile \
-w $outDir \
--min-length $minLength \
$opts all

wait
