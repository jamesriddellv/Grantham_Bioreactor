## Get pileup info ##

# Create output directory
WDIR=/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/02-get-relative-abundance/results/vOTU/read_mapping_metaT/
pileup_dir=${WDIR}/pileup_data
bam_dir=${WDIR}/bam

# Loop through all BAM files
for bam_file in ${bam_dir}/*_SORTED.bam; do
    # Get base name without extension
    base_name=$(basename "$bam_file" _SORTED.bam)

    # Extract reads mapping to contig_591846 and generate pileup
    samtools mpileup -r contig_591846 "$bam_file" > "${pileup_dir}/${base_name}_contig_591846.pileup"
done

## combine into one file ##

# Create a combined file with sample information
echo "sample,position,coverage" > ${pileup_dir}/combined_pileup.csv

for pileup_file in ${pileup_dir}/*.pileup; do
    base_name=$(basename "$pileup_file" _contig_591846.pileup)

    # Extract position and coverage columns, add sample name
    awk -v sample="$base_name" 'BEGIN{OFS=","} {print sample, $2, $4}' "$pileup_file" >> ${pileup_dir}/combined_pileup.csv
done

