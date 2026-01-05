import pandas as pd
import re

# load vOTUs
with open('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined.fna.self-blastn.clusters.tsv', 'r') as f:
    vOTUs = [i.split('\t')[0] for i in f.readlines()]

# load checkv summary file
checkv_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad/quality_summary.tsv', sep='\t')
vOTU_checkv_df = checkv_df.loc[checkv_df['contig_id'].isin(vOTUs)]

# load genomad summary file
genomad_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/all_viral_contigs_summary.tsv', sep='\t')
vOTU_genomad_df = genomad_df.loc[genomad_df['seq_name'].isin(vOTUs)]

# merge the tables
vOTU_checkv_df = vOTU_checkv_df.rename(columns={'contig_id': 'vOTU'})
vOTU_genomad_df = vOTU_genomad_df.rename(columns={'seq_name': 'vOTU'})

vOTU_df = vOTU_genomad_df.merge(vOTU_checkv_df, on='vOTU', how='outer')

# make sure everything is the same length
if (len(vOTU_genomad_df) != len(vOTU_df)) or (len(vOTU_genomad_df) != len(vOTU_checkv_df)):
    raise ValueError(f'genomad, checkv, and merged dataframe are not the same length. Please make sure vOTU headers match.')

# get provirus designation
def is_provirus(row):
    if 'provirus' in row['vOTU']:
        return True
    elif row['provirus'] == 'Yes':
        return True
    else:
        return False

vOTU_df['is_provirus_aggregated'] = vOTU_df.apply(is_provirus, axis=1)

# load dramv
vOTU_dramv_df = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/DRAMv/DRAMv-annotate/annotations.tsv', sep='\t')
vOTU_dramv_df = vOTU_dramv_df.rename(columns={'scaffold': 'vOTU'})
vOTU_dramv_df['vOTU'] = vOTU_dramv_df['vOTU'].apply(lambda x: x.split('-cat')[0])
vOTU_dramv_df['vOTU'] = vOTU_dramv_df['vOTU'].apply(lambda x: x.replace('_provirus', '|provirus'))
vOTU_dramv_df['gene_position'] = vOTU_dramv_df['gene_position'].astype(str)

# there are two vOTUs without dramv annotations, but they also have incredibly high contamination.
# they are present in the DRAMv input file, but not in the scaffolds.fna created by dramv. 
# We will ignore these and not consider them in our final vOTU list.

# merge dramv to merged dataframe
vOTU_df = vOTU_df.merge(vOTU_dramv_df, on='vOTU', how='left')

# genomad annotations
genomad_annotations = pd.read_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/genomad_annotate/combined.fna.self-blastn.clusters_annotate/combined.fna.self-blastn.clusters_genes.tsv', sep='\t')
genomad_annotations = genomad_annotations[['gene','annotation_accessions', 'annotation_description', 'evalue']]
genomad_annotations[['vOTU', 'gene_position']] = genomad_annotations['gene'].str.rsplit('_', n=1, expand=True)
genomad_annotations = genomad_annotations[['vOTU', 'gene_position', 'annotation_accessions', 'annotation_description', 'evalue']]

# merge genomad annotations to merged dataframe
vOTU_df = vOTU_df.merge(genomad_annotations, on=['vOTU', 'gene_position'], how='left')

# fill NAs with ''
vOTU_df = vOTU_df.fillna('')

# get 100 random vOTUs
import random

random.seed(42)

manual_vOTU_check_list = random.sample(vOTUs, 100)
manual_vOTU_check_df = vOTU_df.loc[vOTU_df['vOTU'].isin(manual_vOTU_check_list)]
manual_vOTU_check_df.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTU_manual_check_100.tsv', sep='\t', index=False)

# six vOTUs from manual inspection were marked as "not viral". 
# We will use these to construct rules for not viral sequences.

# use concatenated annotations for filtering

vOTU_df['concat_annotations'] = (
    vOTU_df['kegg_hit'] + ';'
    + vOTU_df['viral_hit'] + ';'
    + vOTU_df['pfam_hits'] + ';'
    + vOTU_df['vogdb_hits'] + ';'
    + vOTU_df['annotation_description']
).fillna('').astype(str)

vOTU_agg = vOTU_df.groupby('vOTU').agg({'concat_annotations': ';'.join}).reset_index()
vOTU_agg = vOTU_agg.merge(vOTU_checkv_df, on='vOTU', how='left')

# clean up with new rules
not_viral_list = []

def manual_filter(vOTU):


    # 1. checkv contamination >80%
    if float(vOTU['contamination']) > 80:
        not_viral_list.append(vOTU['vOTU'])

    # 2. transposase detected, but no checkv viral genes detected
    if 'transposase' in vOTU['concat_annotations'] and int(vOTU['viral_genes']) == 0:
        if vOTU['vOTU'] not in not_viral_list:  # Check by vOTU ID, not the whole dictionary
            not_viral_list.append(vOTU['vOTU'])

    # 3. no checkv viral genes and host genes >30% of genes
    if (str(vOTU['warnings']) == 'no viral genes detected') and (int(vOTU['host_genes']) / int(vOTU['gene_count']) > 0.3):
        if vOTU['vOTU'] not in not_viral_list:
            not_viral_list.append(vOTU['vOTU'])

vOTU_agg.apply(manual_filter, axis=1)

print('number of vOTUs marked as not viral by manual filter: ' + str(len(not_viral_list)))

# extract not viral vOTUs for manual inspection double-check
not_viral_vOTUs = vOTU_df.loc[vOTU_df['vOTU'].isin(not_viral_list)]
cols_to_keep = [
    "vOTU", "topology", "virus_score", "fdr", "n_hallmarks", "contig_length", 
    "provirus", "gene_count", "viral_genes", "host_genes", "checkv_quality", 
    "miuvig_quality", "completeness", "completeness_method", "contamination", 
    "warnings", "is_provirus_aggregated", "gene_position", "strandedness", 
    "kegg_hit", "viral_hit", "pfam_hits", "vogdb_hits", "annotation_description", "is_transposon", 
    "amg_flags"
]

not_viral_vOTUs[cols_to_keep].to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/not_viral_vOTUs_manual_check.tsv', sep='\t', index=False)

# filter and export final vOTU set and metadata

vOTUs_filtered = vOTU_df.loc[~vOTU_df['vOTU'].isin(not_viral_list)]
vOTUs_filtered.to_csv('/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv', sep='\t', index=False)

from Bio import SeqIO

# Input and output file paths
vOTU_fasta_path = "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined.fna.self-blastn.clusters.fna"
output_fasta_path = "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined_manual_filtered.fna"

# List of vOTU headers to keep
vOTUs_filtered_set = set(vOTUs_filtered['vOTU'])  # Convert to a set for faster lookup

# Read and filter the sequences
with open(output_fasta_path, "w") as output_handle:
    for record in SeqIO.parse(vOTU_fasta_path, "fasta"):
        if record.id in vOTUs_filtered_set:
            SeqIO.write(record, output_handle, "fasta")

print(f"Filtered FASTA file saved to {output_fasta_path}")
