#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --mem=180gb
#SBATCH --account=PAS1117
#SBATCH --job-name=recursive-genomad-checkv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=18

######################
### Load variables ###
######################

sample=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ../assemblies.conf)
dbDir=../data/genomad_db

# Define base directories
genomad_base=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad
checkv_base=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/checkv_genomad
recursive_genomad_base=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/genomad/recursive

# Initial input file for Genomad
inFile=../data/grantham_assemblies_5kb/${sample}
iteration=1
all_proviruses=false

####################
### Run pipeline ###
####################

while true; do
    echo "Running Genomad iteration $iteration for sample $sample..."

    if [[ $iteration -eq 1 ]]; then
        genomad_out="${genomad_base}/${sample}"
    else
        genomad_out="${recursive_genomad_base}_iter${iteration}/${sample}"
        inFile="${recursive_genomad_base}_iter$((iteration - 1))/${sample}/proviruses_summary/proviruses_virus.fna"
    fi

    mkdir -p $genomad_out

    genomad end-to-end --cleanup --enable-score-calibration --threads 30 $inFile $genomad_out $dbDir

    # Determine which output file to use for CheckV
    sample_base=$(basename "${sample}" .fa)
    if [[ -d "${genomad_out}/${sample_base}_summary" ]]; then
        # First iteration case (genomad found at least some full viruses)
        checkv_in="${genomad_out}/${sample_base}_summary/${sample_base}_virus.fna"
        all_proviruses=false
    else
        # Subsequent iterations (all were identified as proviruses)
        checkv_in="${genomad_out}/proviruses_summary/proviruses_virus.fna"
        all_proviruses=true
    fi

    checkv_out="${checkv_base}_iter${iteration}/${sample}"
    mkdir -p $checkv_out

    echo "Running CheckV iteration $iteration for sample $sample..."
    time checkv end_to_end $checkv_in $checkv_out -t 30

    # Check if new proviruses exist
    provirus_file="$checkv_out/proviruses.fna"
    if [ ! -s "$provirus_file" ]; then
        echo "No more proviruses detected. Finalizing viral contigs..."
        break
    fi

    # Update input for next iteration
    iteration=$((iteration + 1))
done

###########################
### Merge viral contigs ###
###########################

final_viruses=/fs/scratch/Sullivan_Lab/JamesR/Grantham_Bioreactor/01-build-vOTU-database/results/viral_contigs/${sample}_viruses.fna

cat ${checkv_base}_iter*/${sample}/viruses.fna > $final_viruses

echo "Final viral contigs saved to: $final_viruses"
