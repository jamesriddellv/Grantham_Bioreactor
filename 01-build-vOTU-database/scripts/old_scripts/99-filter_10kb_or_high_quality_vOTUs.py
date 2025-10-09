import pandas as pd
import numpy as np

df = pd.read_csv('../data/viral_id_outputs/checkv/ALL-SAMPLES-quality_summary.tsv', sep='\t')

df = df.loc[df['contig_length'] != 'contig_length']
df['contig_length'] = df['contig_length'].astype(int)

# contig must be >10kb or checkv quality high or Complete
quality_contigs = df.loc[(df['contig_length'] > 10000) | ( (df['checkv_quality'] == 'High-quality') | (df['checkv_quality'] == 'Complete') )]
# load vOTU representative headers
with open('../results/vOTU_clusters/vOTU_cluster_headers.txt', 'r') as f:
    headers = [i.lstrip('>').strip('\n') for i in f.readlines()]

filtered_headers = [i for i in headers if i in list(quality_contigs.contig_id)]

with open('../results/vOTU_clusters/checkv_10kb_or_high_quality_vOTUs.txt', 'w') as f:
    for i in filtered_headers:
        f.write(i + '\n')
