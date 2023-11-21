# this gets the alignment rates for each bowtie2 output with slurm ID 27666949

ls -1 *27666949* | while read line; do tail -2 $line | head -1 | grep -oP '\d+\.\d+'; done | sort -n
