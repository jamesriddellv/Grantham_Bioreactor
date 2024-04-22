module load blast/2.13.0+

makeblastdb \
-in /fs/ess/PAS1117/riddell26/Grantham_Bioreactor/01-build-vOTU-database/data/October2023_Grantham/DRAM_v1.4.4_wCAMPER_2304/scaffolds.fna \
-dbtype nucl \
-out ../01-combined-database/hosts_db
