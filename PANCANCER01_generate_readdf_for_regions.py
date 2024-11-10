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

outdir = "/media/hieunguyen/HNSD_mini/outdir"
PROJECT = "TMD450_TCGA_data_analysis"
thres_hypo = 0.3
thres_hyper = 0.6

path_to_main_output = os.path.join(outdir, PROJECT, output_version)
path_to_03_output = os.path.join(path_to_main_output, "03_output")
path_to_read_data = "/media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions_with_readname_TMDfull"
path_to_01_output = os.path.join(path_to_main_output, "PANCANCER01_output")
os.system(f"mkdir -p {path_to_01_output}")

path_to_save_panel = os.path.join( path_to_main_output, "panel")

cpg450df = pd.read_excel(os.path.join(path_to_save_panel, "TMD450_overlapping_TCGA.xlsx"))
cpg450df = cpg450df[cpg450df['overlapTCGA'] == "yes"]
cpg450df = cpg450df.drop_duplicates(subset=['cpg'])

metadata = pd.read_excel("metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx")

##### generate readdf for all samples, not only training samples
# metadata = metadata[metadata["Set"] == "train"]

metadata.head()
metadata.shape

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
        
all_read_files = [item for item in pathlib.Path(path_to_read_data).glob("*.sorted.csv") if item.name.replace(".sorted.csv", "") in metadata["SampleID"].values]

all_countDMP_results = [item for item in pathlib.Path(path_to_03_output).glob("*/countDMPs.xlsx")]
# remove gastric results, not statistically significant
all_countDMP_results = [item for item in all_countDMP_results if "Gastric" not in str(item)]
print(f"Number of samples in this analysis: {len(all_read_files)}")

testdf = pd.DataFrame()

for item in all_countDMP_results:
    tmpdf = pd.read_excel(item)
    tmpdf["Label"] = str(item).split("/")[-2]
    testdf = pd.concat([testdf, tmpdf], axis = 0)
print("Number of TMD450 regions that have been tested by TCGA data: {}".format(testdf.shape[0]))
if "hyper" not in testdf.columns:
    testdf["hyper"] = 0
if "hypo" not in testdf.columns:
    testdf["hypo"] = 0

testdf["hypo_or_hyper"] = testdf[["hyper", "hypo"]].apply(lambda x: "hyper" if x[0] > x[1] else "hypo", axis = 1)
testdf["region_combiname"] = testdf[["Var1", "hypo_or_hyper"]].apply(lambda x: f"{x[0]}_{x[1]}", axis = 1)

testdf.to_csv(os.path.join(path_to_01_output, "testdf_all_regions_for_pan_cancer.csv"), index = False)

all_samples = []
raw_counts = []
in_read_counts = []
for file in tqdm(all_read_files):
    if os.path.isfile(os.path.join(path_to_01_output, file.name.replace(".sorted.csv", ".read_classification.csv"))) == False:
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
        resdf = tmpdf.groupby(["region", "read_classification"]).seq.count().reset_index().pivot_table(index = "region", 
                                                                                                        columns = "read_classification", 
                                                                                                        values = "seq").reset_index().fillna(0)
        ##### get the region type from TCGA test results
        resdf["region_type"] = resdf["region"].apply(lambda x: ",".join(testdf[testdf.Var1 == x].hypo_or_hyper.unique()))
        if "hyper" not in resdf.columns:
            resdf["hyper"] = 0
        if "hypo" not in resdf.columns:
            resdf["hypo"] = 0
        resdf_orig = resdf.copy()

        resdf = resdf.assign(region_type=resdf["region_type"].str.split(",")).explode("region_type")

        ##### assign candi reads equal to number of hypo or hyper reads, depending on the region type
        resdf["candi_reads"] = resdf[["region_type", "hyper", "hypo"]].apply(lambda x: x[1] if x[0] == "hyper" else x[2], axis = 1)
        
        ##### save the results
        resdf.to_csv(os.path.join(path_to_01_output, "{}.candi_reads.csv".format(file.name.split(".")[0])), index = False)
        tmpdf.to_csv(os.path.join(path_to_01_output, file.name.replace(".sorted.csv", ".read_classification.csv")))
    else:
        print("File {} already exists".format(file.name.replace(".sorted.csv", ".read_classification.csv")))
countdf = pd.DataFrame({"SampleID": all_samples, "raw_count": raw_counts, "in_read_count": in_read_counts})
countdf.to_csv(os.path.join(path_to_01_output, "all_count.csv"))

