# Introduction

# Pipeline

## Preparation

- Run `filter_cov_files.sh` to filter and keep only CpG sites that have more than N reads. 

- Run `generate_regions_from_TMD_regions.ipynb` to generate a table containing all CpG sites in 450 TMD regions. 

- Run `compare_TMD450_and_TCGA_CpG_sites.R` to find overlapping CpG sites between the panel TMD 450 regions and the TCGA methylation array panel. 

- Run `filter_overlapCpG_with_TCGA_from_TMD450.R` to get only CpG sites overlapping between TMD450 and TCGA panel.

- For an input `.BAM` file, run the function `fetch_reads_from_450_regions.py` to collect reads overlapping 450 regions and place them all in a dataframe. The first column is the `read ID`, which we can use to group 2 reads of a pair to form a DNA fragment. To run parallely multiple samples, use the script `run_nextflow.nf`.

## QC CpG sites in 450 regions

- Run `02_generate_counts_for_each_minCov.R` to generate the count-QC for all `cov` files w.r.t different threshold on CpG sequencing depths.

