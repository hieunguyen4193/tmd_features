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

data_version = "TMD_cov"
output_version = "20240910"

outdir = "/media/hieunguyen/GSHD_HN01/outdir"
PROJECT = "TMD450_TCGA_data_analysis"
thres_hypo = 0.3
thres_hyper = 0.6

def main(args):
    input_cancer_class = args.input_cancer_class
    mode = args.mode
    path_to_read_data = args.input
    print("working on input cancer class {}".format(input_cancer_class))

    path_to_main_output = os.path.join(outdir, PROJECT, output_version)
    if mode == "all":
        path_to_03_output = os.path.join(path_to_main_output, "03_output", input_cancer_class)
        path_to_13_output = os.path.join(outdir, PROJECT, output_version, "13_output", input_cancer_class, "thres_hypo_{}_hyper_{}".format(thres_hypo, thres_hyper))
    elif mode == "hypo_only":
        path_to_03_output = os.path.join(path_to_main_output, "03_output_all_hypo", input_cancer_class)
        path_to_13_output = os.path.join(outdir, PROJECT, output_version, "13_output_all_hypo", input_cancer_class, "thres_hypo_{}_hyper_{}".format(thres_hypo, thres_hyper))
    elif mode == "hyper_only":
        path_to_03_output = os.path.join(path_to_main_output, "03_output_all_hyper", input_cancer_class)
        path_to_13_output = os.path.join(outdir, PROJECT, output_version, "13_output_all_hyper", input_cancer_class, "thres_hypo_{}_hyper_{}".format(thres_hypo, thres_hyper))
    
    os.system("mkdir -p {}".format(path_to_13_output))

    path_to_save_panel = os.path.join( path_to_main_output, "panel")

    cpg450df = pd.read_excel(os.path.join(path_to_save_panel, "TMD450_overlapping_TCGA.xlsx"))
    cpg450df = cpg450df[cpg450df['overlapTCGA'] == "yes"]
    cpg450df = cpg450df.drop_duplicates(subset=['cpg'])

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
            
    all_read_files = [item for item in pathlib.Path(path_to_read_data).glob("*.sorted.csv") if item.name.replace(".sorted.csv", "")]
    testdf = pd.read_excel(os.path.join(path_to_03_output, "countDMPs.xlsx"))
    print("Number of TMD450 regions that have been tested by TCGA data: {}".format(testdf.shape[0]))
    if "hyper" not in testdf.columns:
        testdf["hyper"] = 0
    if "hypo" not in testdf.columns:
        testdf["hypo"] = 0
    testdf["hypo_or_hyper"] = testdf[["hyper", "hypo"]].apply(lambda x: "hyper" if x[0] > x[1] else "hypo", axis = 1)

    all_samples = []
    raw_counts = []
    in_read_counts = []
    for file in tqdm(all_read_files):
        if os.path.isfile(os.path.join(path_to_13_output, file.name.replace(".sorted.csv", ".read_classification.csv"))) == False:
            tmpdf = pd.read_csv(file, index_col = [0])
            tmpdf["read_overlap_rate"] = tmpdf[["start", "seq", "region"]].apply(lambda x: check_read_inside_region(x[0], x[1], x[2]), axis = 1)
            raw_count = tmpdf.shape[0]
            in_read_count = tmpdf[tmpdf["read_overlap_rate"] == "in"].shape[0]
            all_samples.append(file.name.replace(".sorted.csv", ""))
            raw_counts.append(raw_count)
            in_read_counts.append(in_read_count)

            ##### keep only reads that are completely inside the region
            tmpdf = tmpdf[tmpdf["read_overlap_rate"] == "in"]

            ##### assign read type: hyper or hypo reads based on the given thresholds
            tmpdf["read_classification"] = tmpdf["alpha"].apply(lambda x: assign_read_type(x, thres_hypo, thres_hyper))
            
            ##### considers only regions that are tested with the TCGA data
            tmpdf["region"] = tmpdf["region"].apply(lambda x: x.replace(":", "_").replace("-", "_"))
            tmpdf = tmpdf[tmpdf["region"].isin(testdf.Var1.unique())]
            ##### count hypo and hyper reads in each region
            resdf = tmpdf.groupby(["region", "read_classification"]).seq.count().reset_index().pivot_table(index = "region", columns = "read_classification", values = "seq").reset_index().fillna(0)

            ##### get the region type from TCGA test results
            resdf["region_type"] = resdf["region"].apply(lambda x: testdf[testdf.Var1 == x].hypo_or_hyper.values[0])
            if "hyper" not in resdf.columns:
                resdf["hyper"] = 0
            if "hypo" not in resdf.columns:
                resdf["hypo"] = 0

            ##### assign candi reads equal to number of hypo or hyper reads, depending on the region type
            resdf["candi_reads"] = resdf[["region_type", "hyper", "hypo"]].apply(lambda x: x[1] if x[0] == "hyper" else x[2], axis = 1)
            
            ##### save the results
            resdf.to_csv(os.path.join(path_to_13_output, "{}.candi_reads.csv".format(file.name.split(".")[0])), index = False)
            tmpdf.to_csv(os.path.join(path_to_13_output, file.name.replace(".sorted.csv", ".read_classification.csv")))
        else:
            print("File {} already exists".format(file.name.replace(".sorted.csv", ".read_classification.csv")))
    countdf = pd.DataFrame({"SampleID": all_samples, "raw_count": raw_counts, "in_read_count": in_read_counts})
    print(countdf.shape)
    countdf.to_csv(os.path.join(path_to_13_output, "all_count.csv"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process read data and classify reads.")
    parser.add_argument("--input", type=str, required=True, help="path to read files")
    parser.add_argument("--input_cancer_class", type=str, required=True, help="input cancer class")
    parser.add_argument("--mode", type=str, required=True, help="choose all or hypo or hyper only")
    
    args = parser.parse_args()
    
    main(args)
    

# for input_cancer_class in Breast CRC Liver Lung;do for mode in all hypo_only hyper_only;do \
#     python 13_generate_readdf_for_regions_validations_ViTruong.py \
#     --input /media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions_validations_Vi-Truong \
#     --input_cancer_class $input_cancer_class --mode $mode;done;done

# for mode in all hypo_only hyper_only;do python 13_generate_readdf_for_regions_validations_ViTruong.py --input /media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions_validations_Vi-Truong --input_cancer_class CRC --mode $mode;done;