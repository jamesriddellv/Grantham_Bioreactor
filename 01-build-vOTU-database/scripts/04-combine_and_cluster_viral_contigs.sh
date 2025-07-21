#!/bin/bash

#########################################
### Specify slurm batch job variables ###
#########################################

#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=vOTU-clustering
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


### 1. Combine checkv viruses and proviruses ###

checkvDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad
tmp=${checkvDir}/tmp.txt
viralContigs=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad/combined.fna
outDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs

cat ${checkvDir}/viruses.fna ${checkvDir}/unmodified_proviruses.fna > ${viralContigs}

### 2. Cluster into vOTUs ###

module use /fs/project/PAS1117/modulefiles
module load singularityImages
module load blast

CheckV-0.8.1-ClusterONLY.sif \
-i $viralContigs \
-t 30 \
-o $outDir \
--min-ani 95 \
--min-qcov 0 \
--min-tcov 85 \
-c \
-f
