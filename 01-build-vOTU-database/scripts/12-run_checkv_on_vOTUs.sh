#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=checkv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


checkvLoc="/users/PAS1117/osu9664/eMicro-Apps/CheckV-0.8.1.sif"
outDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_viral_contigs_stordalen_vOTUs_combined/checkv"
inFile="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vs2_sop_viral_contigs_stordalen_vOTUs_combined/vs2_sop_viral_contigs_stordalen_vOTUs_combined.fa.self-blastn.clusters.fna"

module load singularity
time $checkvLoc end_to_end $inFile $outDir -t 40

cat $outDir/proviruses.fna $outDir/viruses.fna > $outDir/combined.fna
