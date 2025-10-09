#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH -n 30
#SBATCH --account=PAS1117
#SBATCH --job-name=vscore
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

mkdir ../results/v-search

cd ../results/v-search

vOTUs="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered.fna"
active_vOTUs="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/active_vOTU_relative_abundance"

# seqkit grep -f ${active_vOTUs}.txt $vOTUs > ${active_vOTUs}.fna

# get pyrodigal

python3 ../../scripts/parallel-prodigal-gv.py \
-t 30 \
-q \
-i ${active_vOTUs}.fna \
-f gff \
-o ${outDir}/${active_vOTUs}_genes.gff \
-a ${outDir}/${active_vOTUs}_genes.faa \
-d ${outDir}/${active_vOTUs}_genes.fna

hmmsearch -o hmmsearch_out_for_viruses_VS_VOG --tblout hmmsearch_out_for_viruses_VS_VOG.hmmout -E 1e-5 /path/to/VOG_hmm_profiles ${active_vOTUs}_genes.faa

mmseqs createdb ${active_vOTUs}_gene_translations.faa mmseqs_query_seq
mmseqs search mmseqs_query_seq path/to/phrogs_profile_db mmseqs_results_mmseqs ./tmp -e 1e-5
mmseqs createtsv mmseqs_query_seq path/to/phrogs_profile_db mmseqs_results_mmseqs mmseqs_results.tsv

