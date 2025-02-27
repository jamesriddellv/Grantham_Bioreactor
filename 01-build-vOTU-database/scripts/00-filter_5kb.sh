#!/bin/bash
input_dir=../data/grantham_assemblies_1kb
filtered_dir=../data/grantham_assemblies_5kb
mkdir -p "${filtered_dir}"

while read -r assembly; do
  # Full input path
  input_path="${input_dir}/${assembly}"

  # Determine the new filename with replacements
  new_filename=$(echo "${assembly}" | sed -e 's/contigs_1000/contigs_5000/' -e 's/contigs_1kb/contigs_5kb/')

  # Filter sequences >= 5000 bp
  seqkit seq -m 5000 "${input_path}" > "${filtered_dir}/${new_filename}"
done < ../assemblies.conf
