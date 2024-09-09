gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
source(file.path(path.to.main.src, "prepare_promoter_regions.R"))
library(GenomicRanges)
source(file.path(path.to.main.src, "configs.R"))

library(dplyr)
library(tidyverse)
library(minfi)
library(comprehenr)
library(limma)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(mltools)
library(DMRcate)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)

if ("liftOver" %in% installed.packages() == FALSE){
  BiocManager::install("liftOver", update = FALSE)  
}
library(liftOver)

data.version <- "TMD_cov"
output.version <- "20240907"

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/tmd_features"

path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version)
path.to.save.QC.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version, "QC")

cpg450df <- read.csv(file.path(path.to.main.src, "all_cpg_in_450_tmd_regions.csv"))
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(cpg = sprintf("%s_%s", chrom, pos)) %>%
  subset(select = -c(X))

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame()

ann450k <- ann450k %>% rowwise() %>%
  mutate(end = pos + 1)
ann450k <- ann450k[c("chr", "pos", "end", "Name")]
ann450k <- makeGRangesFromDataFrame(df = ann450k, seqnames.field = "chr", start.field = "pos", end.field = "end", keep.extra.columns = TRUE)
chain <- import.chain(file.path(path.to.project.src, "hg19ToHg38.over.chain"))
ann450k.hg38 <- liftOver(ann450k, chain) %>% as.data.frame()  %>% subset(select = c(seqnames, start, Name))
colnames(ann450k.hg38) <- c("chr", "pos", "Name")
ann450k.hg38 <- (ann450k.hg38) %>% rowwise() %>%
  mutate(pos = sprintf("%s_%s", chr, pos))

overlap.tmd450.tcga <- intersect(ann450k.hg38$pos,  cpg450df$cpg)
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(overlapTCGA = ifelse(cpg %in% overlap.tmd450.tcga, "yes", "no"))

cpg450df.overlap <- 

