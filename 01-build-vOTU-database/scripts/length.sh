#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 20
#SBATCH -A PAS1117
#SBATCH -J length

perl /fs/project/PAS1117/Global_AMGs/scripts/lengthcont.pl \
../results/vOTU_clusters/vOTU_clusters_5kb_cutoff.fna \
> ../results/COBRA_contig_ext/vOTU_clusters_5kb_cutoff.length
