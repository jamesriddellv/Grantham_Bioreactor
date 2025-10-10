import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# File Paths
filtered_vOTUs_metadata_file = '/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTUs_filtered_metadata.tsv'
cluster_file = '/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/combined.fna.self-blastn.clusters.tsv'
active_vOTU_path = '/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/02-get-relative-abundance/results/metaT/active_1_gene_per_10kb_vOTUs.txt'
output_venn1 = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/S3-vOTU-origins-no-cutoff.pdf'
output_venn2 = '/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/figures/01-C_1_gene_per_10kb_vOTUs.pdf'
output_origin_table = '/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/vOTUs/vOTU_origins.tsv'

# Contig Ranges
contig_ranges = {
    "grantham": (1, 1166227),
    "sun2024": (1166228, 1171278),
    "emerson2018": (1171279, 1173181)
}

# Load filtered vOTUs
def load_vOTUs(file_path):
    with open(file_path, 'r') as f:
        return set(i.strip().split('\t')[0] for i in f.readlines()[1:])

filtered_vOTUs = load_vOTUs(filtered_vOTUs_metadata_file)

# Load vOTU clusters
vOTUs_df = pd.read_csv(cluster_file, header=None, sep='\t', names=['rep', 'members'])
filtered_vOTUs_df = vOTUs_df[vOTUs_df['rep'].isin(filtered_vOTUs)]

# Generate contig lists with "contig_" prepended
contig_lists = {
    key: set(f"contig_{i}" for i in range(start, end + 1))
    for key, (start, end) in contig_ranges.items()
}

# Function to check presence in clusters
def is_in_cluster(lst, df):
    vOTU_set = set(lst)
    return df['members'].apply(lambda cluster: int(bool(set(cluster.split(',')) & vOTU_set)))

# Apply the function to each category
for category, contigs in contig_lists.items():
    filtered_vOTUs_df[category] = is_in_cluster(contigs, filtered_vOTUs_df)

# Get vOTU sets for Venn diagram
def get_vOTU_sets(df):
    return {category: set(df[df[category] >= 1]['rep']) for category in contig_lists}

vOTU_sets = get_vOTU_sets(filtered_vOTUs_df)

# Load active vOTUs and filter
active_vOTUs = load_vOTUs(active_vOTU_path)
active_filtered_vOTUs_df = filtered_vOTUs_df[filtered_vOTUs_df['rep'].isin(active_vOTUs)]
active_vOTU_sets = get_vOTU_sets(active_filtered_vOTUs_df)

# Function to determine origin
def determine_origin(row):
    origins = []
    if row['grantham'] >= 1:
        origins.append('Grantham')
    if row['sun2024'] >= 1:
        origins.append('Sun 2024')
    if row['emerson2018'] >= 1:
        origins.append('Emerson 2018')
    return ';'.join(origins)

# Create origin table
filtered_vOTUs_df['origin'] = filtered_vOTUs_df.apply(determine_origin, axis=1)
origin_table = filtered_vOTUs_df[['rep', 'origin']].rename(columns={'rep': 'vOTU'})

# Save origin table
origin_table.to_csv(output_origin_table, sep='\t', index=False)

# Function to plot Venn diagram
def plot_venn_diagram(sets, title, output_path):
    colors = ["#0072B2", "#E69F00", "#009E73"]  # Colorblind-friendly colors
    plt.figure(figsize=(6,6), dpi=300)

    venn = venn3([sets['grantham'], sets['sun2024'], sets['emerson2018']], 
                 set_labels=('Grantham', 'Sun 2024', 'Emerson 2018'),
                 set_colors=colors, alpha=0.7)

    # Customize borders and labels
    for patch in venn.patches:
        if patch:
            patch.set_edgecolor("black")
            patch.set_linewidth(1.5)

    for label in venn.subset_labels:
        if label:
            label.set_fontsize(12)
            label.set_fontweight('bold')
            label.set_bbox(dict(facecolor='white', edgecolor='none', boxstyle='round,pad=0.2'))

    # Adjust set label positions slightly to prevent overlap
    if venn.set_labels:
        for label in venn.set_labels:
            label.set_y(label.get_position()[1] - 0.05)

    plt.title(title, fontsize=14, fontweight='bold')
    plt.savefig(output_path, dpi=300)
    plt.close()

# Generate and save Venn diagrams
plot_venn_diagram(vOTU_sets, "Venn Diagram of vOTU Overlaps", output_venn1)
plot_venn_diagram(active_vOTU_sets, "Venn Diagram of Active vOTUs", output_venn2)
