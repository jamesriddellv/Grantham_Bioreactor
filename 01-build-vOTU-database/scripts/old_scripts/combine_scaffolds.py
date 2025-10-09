import pandas as pd
import os

# Define the pattern and range of indices
pattern = "scaffolds.part_{:03d}/final-viral-scored.tsv"
indices = range(1, 11)  # Assuming you have files from 001 to 010

# Create an empty DataFrame to store the concatenated data
concatenated_df = pd.DataFrame()

# Loop through the indices and concatenate the files
for index in indices:
    file_path = pattern.format(index)
    if os.path.exists(file_path):
        df = pd.read_csv(file_path, sep='\t')
        concatenated_df = pd.concat([concatenated_df, df])

# Write the concatenated DataFrame to a new file
concatenated_df.to_csv("MAGs-final-viral-scored.tsv", sep='\t', index=False)
