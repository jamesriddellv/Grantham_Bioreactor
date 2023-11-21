workDir="/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/viral_id_outputs/VS2-SOP-screen"

# merge scaffold outputs into one fasta, store in the permanent data directory

cat ${workDir}/scaffolds.part_*/final-viral-scored.fa > ../data/bioreactor-MAGs-final-viral-scored.fa

# merge new assemblies into one fasta

cat ${workDir}/*.contigs_1kb/final-viral-scored.fa > ../data/2023-10-31-assemblies-final-viral-scored.fa
