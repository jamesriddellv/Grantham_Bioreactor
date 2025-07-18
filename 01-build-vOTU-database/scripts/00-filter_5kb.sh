#!/bin/bash
input_dir=../data/grantham_assemblies_1kb
filtered_dir=../data/grantham_assemblies_5kb
mkdir -p "${filtered_dir}"
# Initialize a global counter (start from 0 or read from a previous file)
i=0
# Check if there are already any files in the filtered directory and get the highest contig number
if [ "$(ls -A ${filtered_dir})" ]; then
  i=$(awk -F 'contig_' 'BEGIN{max=0} /^>contig_/ {if ($2 > max) max=$2} END{print max}' "${filtered_dir}"/*)
fi

while read -r assembly; do
  # Full input path
  input_path="${input_dir}/${assembly}"
  # Determine the new filename with replacements
  new_filename=$(echo "${assembly}" | sed -e 's/contigs_1000/contigs_5000/' -e 's/contigs_1kb/contigs_5kb/')
  # Filter sequences >= 5000 bp
  seqkit seq -g -m 5000 "${input_path}" > "${filtered_dir}/${new_filename}"
  # Rename contigs with globally unique numbers and add 'contig_' prefix
  awk -v i="$i" '/^>/ {print ">contig_" ++i; next} {print}' "${filtered_dir}/${new_filename}" > "${filtered_dir}/tmp.txt"
  mv "${filtered_dir}/tmp.txt" "${filtered_dir}/${new_filename}"
  # Update the global counter to the last used value from the current file
  i=$(awk -F 'contig_' '/^>contig_/ {i=$2} END {print i}' "${filtered_dir}/${new_filename}")
done < ../assemblies_1kb.conf
