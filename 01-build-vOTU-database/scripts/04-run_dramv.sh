arrayid=$(($1 + 1))
assembly=$(sed "${arrayid}q;d" ../assemblies.conf)

cd /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor

module use /fs/project/PAS1117/modulefiles
module load DRAM

# Annotate
fasta_input="./viral_id_outputs/vs2-pass2/${assembly}/for-dramv/final-viral-combined-for-dramv.fa"
affi_input="./viral_id_outputs/vs2-pass2/${assembly}/for-dramv/viral-affi-contigs-for-dramv.tab"

# Variables to pass to DRAMv annotate
opts="--skip_trnascan --threads 28 --min_contig_size 1000"
outDir="./viral_id_outputs/DRAMv/${assembly}/DRAMv-annotate"

mkdir -p ./viral_id_outputs/DRAMv/${assembly}

time DRAM-v.py annotate -i $fasta_input -v $affi_input -o $outDir $opts

# Then summarize
time DRAM-v.py distill -i $outDir/annotations.tsv -o "./viral_id_outputs/DRAMv/${assembly}/DRAMv-distill"
