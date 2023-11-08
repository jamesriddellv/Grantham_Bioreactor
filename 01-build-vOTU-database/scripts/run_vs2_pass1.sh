arrayid=$(($1 + 1))
assembly=$(sed "${arrayid}q;d" ../assemblies.conf)

cd /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor

vs2Loc="/users/PAS1117/osu9664/eMicro-Apps/VirSorter2-2.2.3.sif"
dataDir="./scaffolds_1kb_renamed"
outDir="./viral_id_outputs/vs2-pass1/${assembly}"
minLength=5000
opts="--keep-original-seq --include-groups dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae --min-score 0.5 -j 20"

module load singularity
time $vs2Loc run \
-i ${dataDir}/${assembly}.fa \
-w $outDir \
--min-length $minLength \
$opts all --rerun-incomplete
