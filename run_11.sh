# inputdir="/media/hieunguyen/HNSD_mini/outdir/TMD450_TCGA_data_analysis/20240910/reads_from_450_regions_spikein";
inputdir="/media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions_LOD_samples" # path to LOD samples

cancer_class=$1;
for mode in all hypo_only hyper_only;do python 11_generate_readdf_for_regions_LOD_samples.py --input_cancer_class ${cancer_class} --mode ${mode} --input $inputdir;done