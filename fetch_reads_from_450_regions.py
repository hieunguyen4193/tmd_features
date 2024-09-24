import pysam
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
import argparse

##### example input args
# inputbam = "/Users/hieunguyen/data/bam_files/TMD_bam_bismark.soterd.bam"
# bedfile = "./resources/panel.hg19_liftover_to_hg38.bed"
# outputdir = "./output"

def main():
    parser = argparse.ArgumentParser(description='Generate an image matrix.')
    parser.add_argument('--bam', type=str, required=True, help='Path to the input bam file')
    parser.add_argument('--bed', type=str, required=True, help='path to the input bed file of 450 regions')
    parser.add_argument('--output', type=str, required=False, help='path to save output')
    args = parser.parse_args()

    outputdir = args.output
    bedfile = args.bed
    inputbam = args.bam
    
    os.system("mkdir -p {}".format(outputdir))

    beddf = pd.read_csv(bedfile, sep = "\t", header = None)
    beddf.columns = ["chrom", "start", "end", "gene", "V4", "strand"]
    beddf["region"] = beddf[["chrom", "start", "end"]].apply(lambda x: "{}:{}-{}".format(x[0], x[1], x[2]), axis = 1)

    finaldf = pd.DataFrame()
    for region in tqdm(beddf.region.unique()):
        bamfile_obj = pysam.AlignmentFile(inputbam).fetch(region = region)
        reads = []
        for read in bamfile_obj:
            reads.append(read)
        readdf = pd.DataFrame()
        readdf["name"] = [read.to_dict()["name"] for read in reads]
        readdf["chrom"] = [read.to_dict()["ref_name"] for read in reads]
        readdf["start"] = [read.to_dict()["ref_pos"] for read in reads]
        readdf["cigar"] = [read.to_dict()["cigar"] for read in reads]
        readdf["flen"] = [read.to_dict()["length"] for read in reads]
        readdf["seq"] = [read.to_dict()["seq"] for read in reads]
        readdf["methyl_string"] = [read.to_dict()["tags"][2] for read in reads]
        readdf["XR"] = [read.to_dict()["tags"][3] for read in reads]
        readdf["XG"] = [read.to_dict()["tags"][4] for read in reads]
        readdf["num_CpG"] = readdf["methyl_string"].apply(lambda x: len([item for item in x if item in ["Z", "z"]]))
        readdf["num_methyl"] = readdf["methyl_string"].apply(lambda x: len([item for item in x if item in ["Z"]]))
        readdf["num_unmethyl"] = readdf["methyl_string"].apply(lambda x: len([item for item in x if item in ["z"]]))    
        readdf["alpha"] = readdf[["num_methyl", "num_CpG"]].apply(lambda x: x[0]/x[1], axis = 1)
        readdf["region"] = region
        finaldf = pd.concat([finaldf, readdf], axis = 0)

    finaldf.to_csv(os.path.join(outputdir, inputbam.split("/")[-1].replace(".bam", ".csv")))
    
if __name__ == '__main__':
    main()