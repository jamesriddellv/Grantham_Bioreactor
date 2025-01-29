# combines genomad virus predictions from April2023 assemblies, Stordalen MAGs, and October2023 assemblies for the Grantham Bioreactor

cat ../results/viral_id_outputs/genomad/April2023/ALL-SAMPLES-20230406-STM_0716_E_M_megahit_final.contigs_1000_summary/ALL-SAMPLES-20230406-STM_0716_E_M_megahit_final.contigs_1000_virus.fna \
../results/viral_id_outputs/genomad/April2023/ALL-SAMPLES-20230406-STM_0716_E_M_megahit_final.contigs_1000_find_proviruses/ALL-SAMPLES-20230406-STM_0716_E_M_megahit_final.contigs_1000_provirus.fna \
../results/viral_id_outputs/genomad/October2023/October2023_viruses.fna \
../results/viral_id_outputs/genomad/MAGs/ALL-MAG-scaffolds_summary/ALL-MAG-scaffolds_virus.fna \
../results/viral_id_outputs/genomad/MAGs/ALL-MAG-scaffolds_find_proviruses/ALL-MAG-scaffolds_provirus.fna \
> ../results/viral_id_outputs/genomad/April2023_October2023_ALL-MAG_viruses_combined.fna
