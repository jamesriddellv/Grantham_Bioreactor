#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=checkv-on-genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

######################
### load variables ###
######################

inFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/all_viral_contigs.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad
mkdir -p $outDir

time /fs/ess/PAS1117/modules/singularity/CheckV-1.0.3-updated.sif end_to_end $inFile $outDir -t 30 --remove_tmp

# Handle situations where CheckV identifies additional proviruses

echo "number of additional proviruses detected that need to be re-classified by genomad"
grep ">" $outDir/proviruses.fna | wc -l

python3 03-1-unmodify_provirus_headers.py
