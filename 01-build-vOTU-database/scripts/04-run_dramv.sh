arrayid=$(($1 + 1))
assembly=$(sed "${arrayid}q;d" ../April2023_assemblies.conf)

module use /fs/project/PAS1117/modulefiles
module load DRAM

# For April2023 assemblies

# Annotate
fasta_input="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/viral_id_outputs/April2023_viral_identification/viral_id_outputs/vs2-pass2/${assembly}/for-dramv/final-viral-combined-for-dramv.fa"
affi_input="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/viral_id_outputs/April2023_viral_identification/viral_id_outputs/vs2-pass2/${assembly}/for-dramv/viral-affi-contigs-for-dramv.tab"

# Variables to pass to DRAMv annotate
opts="--skip_trnascan --threads 20 --min_contig_size 1000"
outDir="/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/viral_id_outputs/April2023_viral_identification/viral_id_outputs/DRAMv/${assembly}/DRAMv-annotate"

mkdir -p /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/viral_id_outputs/April2023_viral_identification/viral_id_outputs/DRAMv/${assembly}

time DRAM-v.py annotate -i $fasta_input -v $affi_input -o $outDir $opts

# Then summarize
time DRAM-v.py distill -i $outDir/annotations.tsv -o "/fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/viral_id_outputs/April2023_viral_identification/viral_id_outputs/DRAMv/${assembly}/DRAMv-distill"
