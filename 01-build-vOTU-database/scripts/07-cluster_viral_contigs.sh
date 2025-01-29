#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=cluster_viruses
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=%x_%j.out

### genomad

# same as below, but was double-checking to make sure it wasn't just the conservative post-calibration set (this was the conservative set in question)
# inFile=../results/viral_contigs/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated_5kb.fna
# outDir=../results/vOTU_clusters/genomad_vOTU_clusters

# default

# inFile=../results/viral_contigs/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated_default_5kb.fna
# outDir=../results/vOTU_clusters/genomad_vOTU_default_5kb_clusters

### vs2_sop

# inFile=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/vs2_sop_viral_contigs_5kb_or_complete_combined_deduplicated.fa
# outDir=../results/vs2_sop_vOTU_clusters

### vs2_sop, genomad

# cat ../results/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated_5kb.fna ../results/viral_contigs/vs2_sop_viral_contigs_5kb_or_complete_combined_deduplicated.fa > ../results/viral_contigs/vs2_sop_genomad_viral_contigs_combined.fna
# inFile=../results/viral_contigs/vs2_sop_genomad_viral_contigs_combined.fna
# outDir=../results/vs2_sop_genomad_combined

### vs2_sop, stordalen

# cat ../results/viral_contigs/vs2_sop_viral_contigs_5kb_or_complete_combined_deduplicated.fa ../data/stordalen_vOTUs_symlink.fa > ../results/viral_contigs/vs2_sop_viral_contigs_stordalen_vOTUs_combined.fa
# inFile=../results/viral_contigs/vs2_sop_viral_contigs_stordalen_vOTUs_combined.fa
# outDir=../results/vs2_sop_viral_contigs_stordalen_vOTUs_combined

### genomad, stordalen

# cat ../results/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated_5kb.fna ../data/stordalen_vOTUs_symlink.fa > ../results/viral_contigs/genomad_stordalen_vOTUs_combined.fa
# inFile=../results/viral_contigs/genomad_stordalen_vOTUs_combined.fa
# outDir=../results/genomad_stordalen_vOTUs_combined

### vs2_sop, genomad, stordalen

# cat ../results/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated_5kb.fna ../results/viral_contigs/vs2_sop_viral_contigs_5kb_or_complete_combined_deduplicated.fa ../data/stordalen_vOTUs_symlink.fa > ../results/viral_contigs/vs2_sop_genomad_stordalen_combined.fa
# inFile=../results/viral_contigs/vs2_sop_genomad_stordalen_combined.fa
# outDir=../results/vs2_sop_genomad_stordalen_combined_vOTUs

### vs2_sop, genomad, stordalen, emerson_2018

cat ../results/viral_contigs/genomad/April2023_October2023_ALL-MAG_viruses_combined_deduplicated_default_5kb.fna ../results/viral_contigs/vs2_sop_viral_contigs_5kb_or_complete_combined_deduplicated.fa ../data/stordalen_vOTUs_symlink.fa ../data/emerson_2018_vOTUs.fa > ../results/viral_contigs/vs2_sop_genomad_stordalen_emerson2018_combined.fa
inFile=../results/viral_contigs/vs2_sop_genomad_stordalen_emerson2018_combined.fa
outDir=../results/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs

mkdir -p $outDir

module use /fs/project/PAS1117/modulefiles
module load singularityImages

CheckV-0.8.1-ClusterONLY.sif \
-i $inFile \
-t 40 \
-o $outDir \
--min-ani 95 \
--min-qcov 0 \
--min-tcov 85
