cancer_class=$1;
for mode in hypo_only hyper_only;do python 07_generate_readdf_for_regions.py --input_cancer_class ${cancer_class} --mode ${mode};done