#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH -n 20
#SBATCH --mem=92gb
#SBATCH --account=PAS1117
#SBATCH --job-name=iphop
#SBATCH --mail-user=riddell.26@buckeyemail.osu.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --export=ALL

module use /fs/project/PAS1117/modulefiles
module load iPHoP/1.3.2

export PERL5LIB=/fs/project/PAS1117/modules/iPHoP/1.3.2/lib/perl5/site_perl/5.22.0/:$PERL5LIB

iphop predict \
--fa_file ../../01-build-vOTU-database/results/vOTU_clusters/vOTU_clusters_5kb_cutoff.fna \
--db_dir $DB \
--out_dir ../results/iphop_standard_predictions \
-t 20

