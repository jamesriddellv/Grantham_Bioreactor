arrayid=$(($1 + 1))
assembly=$(sed "${arrayid}q;d" ../assemblies.conf)

cd /fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor

checkvLoc="/users/PAS1117/osu9664/eMicro-Apps/CheckV-0.8.1.sif"
outDir="./viral_id_outputs/checkv/${assembly}"
inFile="./viral_id_outputs/vs2-pass1/${assembly}/final-viral-combined.fa"

module load singularity
time $checkvLoc end_to_end $inFile $outDir -t 40

cat $outDir/proviruses.fna $outDir/viruses.fna > $outDir/combined.fna
