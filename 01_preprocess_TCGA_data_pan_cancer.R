gc()
rm(list = ls())

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)

# downgrading the package matrixStats to 1.1.0 to solve the problem
# Error: useNames = NA is defunct. Instead, specify either useNames = TRUE or useNames = FALSE. 
# install.packages("https://cran.r-project.org/src/contrib/Archive/matrixStats/matrixStats_1.1.0.tar.gz", repos = NULL, type = "source")

tcga.data.version <- "pool"
output.version <- "20240910"

path.to.tcga.download.data <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/TCGA_all_idat", tcga.data.version)

path.to.project.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
path.to.tcga.info <- file.path(path.to.project.src, "TCGA_database_info")

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "TMD450_TCGA_data_analysis"
input.cancer.class <- "pan_cancer"
path.to.main.output <- file.path(outdir, PROJECT, output.version)
if (input.cancer.class == "Colon"){
  path.to.01.output <- file.path(path.to.main.output, "01_output", "CRC")    
} else {
  path.to.01.output <- file.path(path.to.main.output, "01_output", input.cancer.class)
}
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

all.normal.files <- Sys.glob(file.path(path.to.tcga.info, sprintf("%s_idat", input.cancer.class), "*Normal_450K.tsv"))
all.tumor.files <- Sys.glob(file.path(path.to.tcga.info, sprintf("%s_idat", input.cancer.class), "*Tumor_450K.tsv"))

sample.sheet.normal <- data.frame()
sample.sheet.cancer <- data.frame()
for (input.file in all.normal.files){
  tmp.sample.sheet.normal <- read.csv(input.file, sep = "\t")  
  sample.sheet.normal <- rbind(sample.sheet.normal, tmp.sample.sheet.normal)
}

for (input.file in all.tumor.files){
  tmp.sample.sheet.cancer <- read.csv(input.file, sep = "\t") 
  sample.sheet.cancer <- rbind(sample.sheet.cancer, tmp.sample.sheet.cancer)
}

sample.sheet <- rbind(sample.sheet.normal, sample.sheet.cancer)

if (nrow(sample.sheet) == 0){
  print(input.cancer.class)
}
print(unique(sample.sheet$Sample.Type))
path.to.pool.data <- file.path(path.to.tcga.download.data, input.cancer.class)
all.idat.files <- Sys.glob(path.to.pool.data)

sample.sheet <- sample.sheet %>% rowwise() %>%
  mutate(path = to_vec(for (item in all.idat.files) if (grepl(File.Name, item)) item)[[1]])
idat.metadata <- data.frame(File.Name = to_vec(
  for (item in sample.sheet$File.Name) str_replace(str_replace(item, "_Red.idat", ""), "_Grn.idat", "")
))

if (file.exists(file.path(path.to.01.output, "idat.obj.rds")) == FALSE){
  idat.obj <- read.metharray.exp(targets = data.frame(Basename = unique(idat.metadata$File.Name)), 
                                 verbose = TRUE, 
                                 force = TRUE, 
                                 base = path.to.pool.data)
  
  saveRDS(idat.obj, file.path(path.to.01.output, "idat.obj.rds"))
} else {
  print("reading in idat object...")
  idat.obj <- readRDS(file.path(path.to.01.output, "idat.obj.rds"))
}

##### DETECTION P VALUE
if (file.exists(file.path(path.to.01.output, "detP.rds")) == FALSE){
  detP <- detectionP(idat.obj)
  saveRDS(detP, file.path(path.to.01.output, "detP.rds"))
} else {
  print("reading in detP ...")
  detP <- readRDS(file.path(path.to.01.output, "detP.rds"))
}

if (file.exists(file.path(path.to.01.output, "idat.obj.preprocessQuantile.rds")) == FALSE){
  # filter high detection probes
  print("filter high detection probes")
  keep <- colMeans(detP) < 0.05
  idat.obj <- idat.obj[,keep]
  
  detP <- detP[,keep]
  
  # preprocessing quantile for idat object
  mSetSq <- preprocessQuantile(idat.obj)  
  
  detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
  
  # remove any probes that have failed in one or more samples
  keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
  mSetSqFlt <- mSetSq[keep,]
  
  # get 450k CpG annotation    
  ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  # remove CpGs/probe in sex chromosomes
  keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
  mSetSqFlt <- mSetSqFlt[keep,]
  
  # remove probes with SNPs at CpG site
  mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
  
  # exclude cross reactive probes 
  xReactiveProbes <- read.csv(file.path(path.to.project.src, "48639-non-specific-probes-Illumina450k.csv") ,sep=",", stringsAsFactors=FALSE)
  keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID)
  mSetSqFlt <- mSetSqFlt[keep,] 
  
  # save object mSetSqFlt
  saveRDS(mSetSqFlt, file.path(path.to.01.output, "mSetSqFlt"))
  print("finished saving mSetSqFlt")
} else {
  print("reading preprocessed quantile data, mSetSq")
  mSetSqFlt <- readRDS(file.path(path.to.01.output, "mSetSqFlt"))
}

mVals <- getM(mSetSqFlt)
bVals <- getBeta(mSetSqFlt)

write.csv(bVals, file.path(path.to.01.output, "bVals.csv"))
saveRDS(bVals, file.path(path.to.01.output, "bVals.rds"))

write.csv(mVals, file.path(path.to.01.output, "mVals.csv"))
saveRDS(mVals, file.path(path.to.01.output, "mVals.rds"))