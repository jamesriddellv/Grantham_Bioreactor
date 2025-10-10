#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -n 10
#SBATCH --account=PAS1117
#SBATCH --job-name=htseq-count
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

outFile="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/htseq_vOTUs_100M_90FILTERED_REVSTRANDED.tsv"
annoFile="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/prodigal-gv/combined_manual_filtered.gff"
bam_dir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/bam"

# run htseq-count
htseq-count -a 0 \
-t CDS \
-i ID \
-n 10 \
--stranded=reverse \
-c "$outFile" \
${bam_dir}/*90FILTERED_NAMESORTED.bam \
"$annoFile"

