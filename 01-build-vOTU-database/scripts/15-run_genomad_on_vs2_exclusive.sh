#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --mem=117gb
#SBATCH --account=PAS1117
#SBATCH --job-name=genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

inFile=../results/vOTU_clusters/vs2_sop_genomad_stordalen_emerson2018_combined_vOTUs/vs2_exclusive_vOTUs.fna
# outDir=../results/viral_id_outputs/genomad
outDir=../results/viral_id_outputs/genomad/vs2_exclusive
dbDir=../data/genomad_db
mkdir -p $outDir

genomad end-to-end --cleanup --enable-score-calibration --threads 28 $inFile $outDir $dbDir
# genomad annotate --restart --cleanup --threads 40 $inFile $outDir $dbDir
# genomad find-proviruses --restart --cleanup --threads 40 $inFile $outDir $dbDir
# genomad marker-classification --restart --cleanup --threads 40 $inFile $outDir $dbDir
# genomad nn-classification --restart --cleanup --threads 40 $inFile $outDir
# genomad score-calibration --restart --cleanup --threads 40 $inFile $outDir
# genomad summary --restart --cleanup --threads 40 $inFile $outDir
