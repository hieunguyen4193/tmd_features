import pandas as pd
import numpy as np
import pathlib 
import os
import matplotlib.pyplot as plt
import seaborn as sns
import random
from tqdm import tqdm
import warnings
import pandas as pd
import argparse
warnings.filterwarnings('ignore')

thres_hypo = 0.3
thres_hyper = 0.6

def assign_read_type(x, thres_hypo, thres_hyper):
    if x < thres_hypo:
        return "hypo"
    elif x > thres_hyper:
        return "hyper"
    else:
        return "none"
def check_read_inside_region(start, seq, region):
        read_end = start + len(seq)
        region_start = int(region.split(":")[1].split("-")[0])
        region_end = int(region.split(":")[1].split("-")[1])
        if start >= region_start and read_end <= region_end:
            return "in"
        else: 
            return "overlap"
        
def main(args):
    path_to_testdf = args.testdf
    path_to_read_data =  args.input
    testdf = pd.read_csv(path_to_testdf)
    outputdir = args.output
    os.system(f"mkdir -p {outputdir}")
    filename = path_to_read_data.split("/")[-1].replace(".sorted.csv", ".read_classification.csv")
    if os.path.isfile(os.path.join(outputdir, filename)) == False:
        tmpdf = pd.read_csv(path_to_read_data, index_col = [0])
        tmpdf["read_overlap_rate"] = tmpdf[["start", "seq", "region"]].apply(lambda x: check_read_inside_region(x[0], x[1], x[2]), axis = 1)
        
        ##### keep only reads that are completely inside the region
        tmpdf = tmpdf[tmpdf["read_overlap_rate"] == "in"]

        ##### assign read type: hyper or hypo reads based on the given thresholds
        tmpdf["read_classification"] = tmpdf["alpha"].apply(lambda x: assign_read_type(x, thres_hypo, thres_hyper))

        ##### considers only regions that are tested with the TCGA data
        tmpdf["region"] = tmpdf["region"].apply(lambda x: x.replace(":", "_").replace("-", "_"))
        tmpdf = tmpdf[tmpdf["region"].isin(testdf.Var1.unique())]
        ##### count hypo and hyper reads in each region
        resdf = tmpdf.groupby(["region", "read_classification"]).seq.count().reset_index().pivot_table(index = "region", 
                                                                                                        columns = "read_classification", 
                                                                                                        values = "seq").reset_index().fillna(0)
        ##### get the region type from TCGA test results
        resdf["region_type"] = resdf["region"].apply(lambda x: ",".join(testdf[testdf.Var1 == x].hypo_or_hyper.unique()))
        if "hyper" not in resdf.columns:
            resdf["hyper"] = 0
        if "hypo" not in resdf.columns:
            resdf["hypo"] = 0
        
        resdf = resdf.assign(region_type=resdf["region_type"].str.split(",")).explode("region_type")

        ##### assign candi reads equal to number of hypo or hyper reads, depending on the region type
        resdf["candi_reads"] = resdf[["region_type", "hyper", "hypo"]].apply(lambda x: x[1] if x[0] == "hyper" else x[2], axis = 1)
        
        ##### save the results
        resdf.to_csv(os.path.join(outputdir, "{}.candi_reads.csv".format(filename.split(".")[0])), index = False)
        tmpdf.to_csv(os.path.join(outputdir, filename.replace(".sorted.csv", ".read_classification.csv")))
    else:
        print(f"File {os.path.join(outputdir, filename)} exists")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process read data and classify reads.")
    parser.add_argument("--input", type=str, required=True, help="input read file")
    parser.add_argument("--output", type=str, required=True, help="output read file")
    parser.add_argument("--testdf", type=str, required=True, help="input testdf file containing all regions hyper and hpo")
    
    args = parser.parse_args()
    
    main(args)


##### run for all TMD samples
# parallel -j 50 python PANCANCER01_generate_readdf_in_parallel.py \
# --input {} \
# --testdf /media/hieunguyen/HNSD_mini/outdir/TMD450_TCGA_data_analysis/20240910/PANCANCER01_output/testdf_all_regions_for_pan_cancer.csv \
# --output /media/hieunguyen/HNSD_mini/outdir/TMD450_TCGA_data_analysis/20240910/PANCANCER01_output ::: $files

##### run for all LOD samples
# path_to_input="/media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions_LOD_samples"
# parallel -j 50 python PANCANCER01_generate_readdf_in_parallel.py \
# --input {} \
# --testdf /media/hieunguyen/HNSD_mini/outdir/TMD450_TCGA_data_analysis/20240910/PANCANCER01_output/testdf_all_regions_for_pan_cancer.csv \
# --output /media/hieunguyen/HNSD_mini/outdir/TMD450_TCGA_data_analysis/20240910/PANCANCER03_output ::: $files

##### run for all validation samples
# path_to_input="/media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions_validations_Vi-Truong"
# files=$(ls ${path_to_input}/*.csv);
# parallel -j 50 python PANCANCER01_generate_readdf_in_parallel.py \
# --input {} \
# --testdf /media/hieunguyen/GSHD_HN01/outdir/TMD450_TCGA_data_analysis/20240910/PANCANCER01_output/testdf_all_regions_for_pan_cancer.csv \
# --output /media/hieunguyen/GSHD_HN01/outdir/TMD450_TCGA_data_analysis/20240910/PANCANCER05_output ::: $files