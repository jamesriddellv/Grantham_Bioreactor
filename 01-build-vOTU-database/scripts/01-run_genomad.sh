#!/bin/bash

### Specify slurm batch job variables ###

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=177gb
#SBATCH --account=PAS1117
#SBATCH --job-name=genomad
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%j.out

######################
### load variables ###
######################


# inFile=../data/ALL-SAMPLES-20230406-STM_0716_E_M_megahit_final.contigs_1000.fa renamed to April2023_contigs_1kb.fa
inFile=../data/October2023_contigs_1kb.fa
# outDir=../results/viral_id_outputs/genomad
outDir=../results/viral_contigs/genomad/October2023
dbDir=../data/genomad_db
mkdir -p $outDir

###################
### run genomad ###
###################

# genomad end-to-end --cleanup --enable-score-calibration --threads 40 $inFile $outDir $dbDir
# genomad annotate --restart --cleanup --threads 40 $inFile $outDir $dbDir
genomad find-proviruses --restart --cleanup --threads 40 $inFile $outDir $dbDir
genomad marker-classification --restart --cleanup --threads 40 $inFile $outDir $dbDir
genomad nn-classification --restart --cleanup --threads 40 $inFile $outDir
genomad score-calibration --restart --cleanup --threads 40 $inFile $outDir
genomad summary --restart --cleanup --threads 40 $inFile $outDir
