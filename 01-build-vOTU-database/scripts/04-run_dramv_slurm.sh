#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH -n 20
#SBATCH --mem=92gb
#SBATCH --account=PAS1117
#SBATCH --job-name=dramv
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --array=1-29%10

bash ./04-run_dramv.sh ${SLURM_ARRAY_TASK_ID}
