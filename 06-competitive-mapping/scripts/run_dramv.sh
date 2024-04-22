module use /fs/project/PAS1117/modulefiles
module load DRAM

# Annotate
fasta_input="../01-combined-database/filtered_clusters.fna"
affi_input="../01-combined-database/for-dramv/viral-affi-contigs-for-dramv.tab"

# Variables to pass to DRAMv annotate
opts="--skip_trnascan --threads 28 --min_contig_size 1000"
outDir=".,/01-combined-database/DRAMv/DRAMv-annotate"

mkdir -p ./01-combined-database/DRAMv/

time DRAM-v.py annotate -i $fasta_input -v $affi_input -o $outDir $opts

# Then summarize
time DRAM-v.py distill -i $outDir/annotations.tsv -o "../01-combined-database/DRAMv/DRAMv-distill"
