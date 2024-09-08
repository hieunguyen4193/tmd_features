import pandas as pd
import numpy as np
import pathlib 
import os
import matplotlib.pyplot as plt
import seaborn as sns

from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

path_to_main_src = pathlib.Path("/media/hieunguyen/HNSD01/src/tmd_features")
data_version = "TMD_cov"
output_version = "20240907"

outdir = pathlib.Path("/media/hieunguyen/HNSD_mini/outdir")
path_to_input = outdir / "raw_data" / "bismark_cov" / data_version
path_to_main_output = outdir / "TMD_read_based_features" / "output" / f"data_{data_version}" / output_version
path_to_difftest_output = path_to_main_output / "difftest_output"
path_to_difftest_output.mkdir(parents=True, exist_ok=True)

path_to_save_QC_output = path_to_main_output / "QC"

thres_hypo = 0.3
thres_hyper = 0.6

input_cancer_class = "Liver"
check_region = pd.read_excel(os.path.join(path_to_difftest_output, "diff_test_{}_region_summary.xlsx".format(input_cancer_class)))

metadata = pd.read_excel("metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx")

path_to_read_data = "/media/hieunguyen/GSHD_HN01/raw_data/reads_from_450_regions"

path_to_py01_output = os.path.join(path_to_main_output, "py01_output")
os.system("mkdir -p {}".format(path_to_py01_output))

all_read_files = [item for item in pathlib.Path(path_to_read_data).glob("*.csv")]

def classify_read(alpha, region):
    region_type = check_region[check_region["CpG"] == region]["region_type"].values[0]
    if region_type == "hypo":
        if alpha < thres_hypo:
            return "candi"
        else:
            return "none"
    elif region_type == "hyper":
        if alpha > thres_hyper:
            return "candi"
        else:
            return "none"

def assign_read_type(x, thres_hypo, thres_hyper):
    if x < thres_hypo:
        return "hypo"
    elif x > thres_hyper:
        return "hyper"
    else:
        return "none"

for file in tqdm(all_read_files):
    tmpdf = pd.read_csv(file, sep=",", index_col=0)

    tmpdf["region"] = tmpdf["region"].apply(lambda x: x.replace(":", "_").replace("-", "_"))
    tmpdf = tmpdf[tmpdf["region"].isin(check_region.CpG.unique())]

    tmpdf["read_type"] = tmpdf["alpha"].apply(lambda x: assign_read_type(x, thres_hypo, thres_hyper))

    tmpdf["read_classification"] = tmpdf[["alpha", "region"]].apply(lambda x: classify_read(x[0], x[1]), axis = 1)

    count_candi_reads = tmpdf.groupby(["region", "read_classification"])["methyl_string"].count().reset_index().pivot(index = "region", columns = "read_classification", values = "methyl_string").fillna(0)
    count_candi_reads["fraction"] = count_candi_reads[["candi", "none"]].apply(lambda x: x[0]/(x[0] + x[1]), axis = 1)

    count_candi_reads.to_csv(os.path.join(path_to_py01_output, "{}.candi_reads.csv".format(file.name.split(".")[0])))