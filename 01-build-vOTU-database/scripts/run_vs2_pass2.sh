arrayid=$(($1 + 1))
assembly=$(sed "${arrayid}q;d" ../assemblies.conf)

cd /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
outDir="./viral_id_outputs/vs2-pass2/${assembly}"
minLength=5000
opts="--seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -j 20 --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-score 0.5"

module load singularity
time $vs2Loc run \
-i ./viral_id_outputs/checkv/${assembly}/combined.fna \
-w $outDir \
--min-length $minLength \
$opts all

wait
