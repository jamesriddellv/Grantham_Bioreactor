import pandas as pd
import numpy as np

ANI="97"
mode="mean"
ALIGNED_PERC="75"
MIN_FRAC="75"

import os
data_dir = f'/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/MAGs/metaG/coverm'  # Replace with the path to your directory containing the files
merged_df = []
for filename in os.listdir(data_dir):
    if filename.endswith('.txt'):  # Assuming the files are in CSV format
        file_path = os.path.join(data_dir, filename)
        df = pd.read_csv(file_path, sep='\t')
        if type(merged_df) == list:
            merged_df = df
        else:
            merged_df = merged_df.merge(df, on='Genome', how='outer')

merged_df.to_csv(data_dir + f'/metaG_MAG_{mode}_{ANI}_{ALIGNED_PERC}_{MIN_FRAC}.tsv', index=False, sep='\t')
