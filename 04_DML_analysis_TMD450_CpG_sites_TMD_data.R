gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)

if ("liftOver" %in% installed.packages() == FALSE){
  BiocManager::install("liftOver", update = FALSE)
}
library(liftOver)

output.version <- "20240910"
data.version <- "TMD_cov"
# install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz", type = "source", repos = NULL)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

path.to.project.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))

path.to.tcga.info <- file.path(path.to.project.src, "TCGA_database_info")

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "TMD450_TCGA_data_analysis"
input.cancer.class <- "Liver"
min.cov <- 5

path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, PROJECT, output.version)
path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("%s_vs_Control_minCov_%s", input.cancer.class, min.cov))
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

all.cov.files <- Sys.glob(file.path(path.to.input, sprintf("filtered_%sreads_TMD450regionsOverlapTCGA_cov", min.cov), "*.cov"))
names(all.cov.files) <- unlist(lapply(all.cov.files, function(x){
  x <- basename(x)
  x <- str_split(x, ".deduplicated")[[1]][[1]]
  x <- str_replace(x, ".cov", "")
  return(x)
}))

##### filter and keep control and the specified cancer class in training class only.
meta.data <- readxl::read_excel(file.path(path.to.project.src, "metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx"))
meta.data <- subset(meta.data, meta.data$Label %in% c("Control", input.cancer.class))
meta.data <- meta.data %>% rowwise() %>%
  mutate(cov.name = str_split(basename(TM_COV), ".deduplicated")[[1]][[1]])
meta.data <- subset(meta.data, is.na(meta.data$cov.name) == FALSE)
meta.data <- subset(meta.data, meta.data$Set == "train")
meta.data <- meta.data[!duplicated(meta.data$cov.name), ]
all.cov.files <- all.cov.files[meta.data$cov.name]
labels <- to_vec( for(item in names(all.cov.files)) if (subset(meta.data, meta.data$cov.name == item)$Label == input.cancer.class) 1 else 0)

##### RUN DML

if (file.exists(file.path(path.to.04.output, "DML_obj.rds")) == FALSE){
  print("Generate methylkit object from input cov files...")
  DML.obj <- readBismarkCoverage( all.cov.files,
                                  sample.id = names(all.cov.files),
                                  assembly = "hg38",
                                  treatment = labels,
                                  context = "CpG",
                                  min.cov = min.cov)
  saveRDS(DML.obj, file.path(path.to.04.output, "DML_obj.rds"))  
} else {
  print("Methylkit object exists, reading in...")
  DML.obj <-  readRDS(file.path(path.to.04.output, "DML_obj.rds"))
}

#####----------------------------------------------------------------------#####
# Unite the meth objects
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.04.output, "meth.rds")) == FALSE){
  print("Unite the methylation methylkit object...")
  meth <- methylKit::unite(object = DML.obj, destrand = FALSE)
  saveRDS(meth, file.path(path.to.04.output, "meth.rds"))  
} else {
  print("Methylkit object united, reading in...")
  meth <- readRDS(file.path(path.to.04.output, "meth.rds"))
}

#####----------------------------------------------------------------------#####
# Generate PCA from CpG loci
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.04.output, "pca_from_CpG_loci_res.rds")) == FALSE){
  print("Generate PCA from CpG methylation loci information...")
  pca.res <- PCASamples(meth, obj.return = TRUE, screeplot = FALSE)
} else {
  print("PCA done, reading in...")
  pca.res <- readRDS(file.path(path.to.04.output, "pca_from_CpG_loci_res.rds"))
}

#####----------------------------------------------------------------------#####
# Generate DML
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.04.output, "DML.rds")) == FALSE){
  print("Conduct differential methylated test between case and control ...")
  myDiff <- calculateDiffMeth(meth, mc.cores = 45)
  saveRDS(myDiff, file.path(path.to.04.output, "DML.rds"))
} else {
  print("Differential methylated test results existed, reading in...")
  myDiff <- readRDS(file.path(path.to.04.output, "DML.rds"))
}

# and get all diff loci
if (file.exists(file.path(path.to.04.output, "diff_locidf.rds")) == FALSE){
  diff.loci <- getMethylDiff(myDiff, difference = methdiff.cutoff, qvalue = qvalue.cutoff)
  saveRDS(diff.loci, file.path(path.to.04.output, "diff_locidf.rds"))  
} else {
  diff.loci <- readRDS(file.path(path.to.04.output, "diff_locidf.rds"))
}
