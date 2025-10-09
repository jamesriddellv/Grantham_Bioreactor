#!/bin/bash
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=prep-genomad-for-dramv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out


vOTU_fasta=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined.fna.self-blastn.clusters.fna
wDir=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/mvip
metadata=../data/metadata_for_mvip.txt
proteins=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/genomad_annotate/combined.fna.self-blastn.clusters_annotate/combined.fna.self-blastn.clusters_proteins.faa
genomad_annotations=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/genomad_annotate/combined.fna.self-blastn.clusters_annotate/combined.fna.self-blastn.clusters_genes.tsv

### make metadata file for MVP pipeline ###

awk 'BEGIN {OFS="\t"; print "Sample_number", "Sample", "Assembly_Path", "Read_Path"} {split($0, a, "/"); print NR, a[length(a)], "../data/grantham_assemblies_5kb/" $0, "../data/grantham_assemblies_5kb/" $0}' ../assemblies_5kb.conf > ${metadata}
mkdir -p ${wDir}

### Run MVP pipeline ###

# mamba activate mvip

# run MVP setup

mvip MVP_00_set_up_MVP \
-i ${wDir} \
-m ${metadata} \
--skip_check_errors


# trick module 06 by generating input files from our outputs without running the whole pipeline

mkdir -p ${wDir}/03_CLUSTERING
cp ${vOTU_fasta} ${wDir}/03_CLUSTERING/MVP_03_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Sequences.fna

mkdir -p ${wDir}/06_FUNCTIONAL_ANNOTATION
cp ${proteins} ${wDir}/06_FUNCTIONAL_ANNOTATION/MVP_06_All_Sample_Filtered_Relaxed_Representative_Virus_Provirus_Protein_Sequences.faa
cp ${genomad_annotations} ${wDir}/06_FUNCTIONAL_ANNOTATION/MVP_06_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Gene_Annotation_GENOMAD.tsv

# rename headers in MVP_06_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Gene_Annotation_GENOMAD.tsv

sed -i '1s/gene/Viral_gene_ID/; 1s/marker/GENOMAD_marker/; 1s/evalue/GENOMAD_evalue/; 1s/bitscore/GENOMAD_Score/; 1s/plasmid_hallmark/GENOMAD_plasmid_hallmark/; 1s/virus_hallmark/GENOMAD_virus_hallmark/' ${wDir}/06_FUNCTIONAL_ANNOTATION/MVP_06_All_Sample_Filtered_Relaxed_Merged_Genomad_CheckV_Representative_Virus_Proviruses_Gene_Annotation_GENOMAD.tsv

mvip MVP_06_do_functional_annotation -i ${wDir} -m ${metadata} --DRAM --threads 30
