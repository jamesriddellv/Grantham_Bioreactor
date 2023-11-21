#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=cluster_viruses
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=%x_%j.out

# cat ../data/*.fa > /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/viral_id_outputs/all_seqs.fna

module use /fs/project/PAS1117/modulefiles
module load singularityImages

CheckV-0.8.1-ClusterONLY.sif \
-i /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/viral_id_outputs/all_seqs.fna \
-t 28 \
-o ../results/vOTU_clusters \
--min-ani 95 \
--min-qcov 0 \
--min-tcov 85
