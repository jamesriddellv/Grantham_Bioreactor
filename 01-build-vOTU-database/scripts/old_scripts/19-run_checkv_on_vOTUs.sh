#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=checkv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

outDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/checkv"
inFile="../results/vOTU_clusters/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vs2_sop_genomad_stordalen_emerson2018_combined.fa.self-blastn.clusters.fna"

time checkv end_to_end $inFile $outDir -t 28

cat $outDir/proviruses.fna $outDir/viruses.fna > $outDir/combined.fna
