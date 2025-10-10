import pandas as pd
import numpy as np

ANI="90"
mode="trimmed_mean"
ALIGNED_PERC="75"
MIN_FRAC="10"

import os
data_dir = f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaG/coverm_trimmed_mean'
merged_df = []
for filename in os.listdir(data_dir):
    if filename.startswith('STM'):  # Assuming the files are in CSV format
        file_path = os.path.join(data_dir, filename)
        df = pd.read_csv(file_path, sep='\t')
        if type(merged_df) == list:
            merged_df = df
        else:
            merged_df = merged_df.merge(df, on='Contig', how='outer')

merged_df.to_csv(data_dir + f'/metaG_vOTU_{mode}_{ANI}_{ALIGNED_PERC}_{MIN_FRAC}.tsv', index=False, sep='\t')
