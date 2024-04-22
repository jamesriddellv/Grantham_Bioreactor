vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
outDir="../01-combined-database/vs2-pass2"
minLength=5000
opts="--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -j 28 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-score 0.5"

module load singularity
time $vs2Loc run \
-i ../01-combined-database/filtered_clusters.fna \
-w $outDir \
--min-length $minLength \
$opts all

wait
