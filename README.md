# Grantham Bioreactor Virus Story

## Purpose

This experiment’s goal is to determine the viral diversity of a bioreactor, how it compares to the original sample site (where the core to generate the bioreactor came from), and how the viral diversity changes in response to perturbation by catechin, a tannin derivative that reduces methanogenesis by competing for available hydrogen (the energy source for methanogenesis)

**Driving Questions**

1. What is the viral diversity of a bioreactor? (alpha diversity)
    1. Which vOTUs are most dominant?
    2. Does the viral diversity change over time during steady state? Are there any clear patterns?
2. How does the viral diversity compare to the Stordalen fen? (beta diversity)
    1. Which key organisms are different? Are the same vOTUs dominant?
3. How does the perturbation by catechin change the composition of the viral community, if at all? Are these viruses at all linked to hosts also impacted by catechin?
4. (potentially own paper / story) Can we metabolically model the bioreactor?

****************Approach****************

1. Build Grantham Bioreactor vOTU database
    1. Mine new assemblies for viruses
    2. Combine with previous bioreactor assemblies
    3. Combine all bioreactor assemblies with Stordalen viruses
    4. Cluster at 95% ANI 80% cov
    5. Use COBRA to extend all vOTUs (maybe?)
        1. Test a subset of vOTUs and the sequences in their cluster and see if when we extend, some of the non-representative sequences extend longer than the vOTUs (which are theoretically the longest).
2. Get relative abundance of vOTUs through time
    1. Map each metaT sample to the vOTU database (bowtie2)
    2. Compute relative abundance (CoverM contig)
    3. Compare vOTUs through time (find papers that do this!)
    4. Plot perturbation line and see if there are any distinct changes
    5. Can also compute beta diversity of pre- and post- perturbation
3. See if any viruses have a significant impact on methane production
    1. Build a WGCNA for methane production and see if any modules are strongly correlated, then pull out the VIPs.

**************************************Directory Structure**************************************

Grantham_Bioreactor

- 01-build-vOTU-database
    - data
        - Stordalen-vOTUs.fna
        - processed-metaG-assemblies
        - new-metaG-assemblies
            - 16 assemblies
        - MAG-database
    - results
        - bioreactor-vOTU-database.fna
- results
- 02-get-relative-abundance
    - data
        - MetaT-reads
            - 1.2Tb metaT data
    - results
        - bam-files
        - CoverM
        - relative-abundance-table
        - WGCNA
    - figures
        - relative-abundance.png
        - alpha-diversity.png
        - beta-diversity-perturbation.png
        - WGCNA.png
- 03-compare-to-Stordalen
    - data
        - Stordalen-relative-abundance
    - results
        - diversity
