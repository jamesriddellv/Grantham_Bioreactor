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
outDir="../results/viral_id_outputs/genomad/checkv/"
inFile="../results/viral_id_outputs/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated.fna"

module load singularity
time $checkvLoc end_to_end $inFile $outDir -t 40

# cat $outDir/proviruses.fna $outDir/viruses.fna > $outDir/combined.fna
