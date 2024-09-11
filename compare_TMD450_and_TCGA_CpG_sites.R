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

# install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/IlluminaHumanMethylation450kanno.ilmn12.hg19_0.6.1.tar.gz", type = "source", repos = NULL)
library("IlluminaHumanMethylation450kanno.ilmn12.hg19")

path.to.project.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
path.to.tcga.info <- file.path(path.to.project.src, "TCGA_database_info")

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "TMD450_TCGA_data_analysis"

path.to.main.output <- file.path(outdir, PROJECT, output.version)
path.to.project.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
path.to.save.panel.info <- file.path(path.to.main.output, "panel")
dir.create(path.to.save.panel.info, showWarnings = FALSE, recursive = TRUE)

cpg450df <- read.csv(file.path(path.to.project.src, "all_cpg_in_450_tmd_regions.csv"))
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(cpg = sprintf("%s_%s", chrom, pos)) %>%
  subset(select = -c(X))

ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) %>% as.data.frame()

ann450k <- ann450k %>% rowwise() %>%
  mutate(end = pos + 1)
ann450k <- ann450k[c("chr", "pos", "end", "Name")]
ann450k <- makeGRangesFromDataFrame(df = ann450k, seqnames.field = "chr", start.field = "pos", end.field = "end", keep.extra.columns = TRUE)

##### lift over the ann450k in hg19 to hg38
chain <- import.chain(file.path(path.to.project.src, "hg19ToHg38.over.chain"))
ann450k.hg38 <- liftOver(ann450k, chain) %>% as.data.frame()  %>% subset(select = c(seqnames, start, Name))
colnames(ann450k.hg38) <- c("chr", "pos", "Name")
ann450k.hg38 <- (ann450k.hg38) %>% rowwise() %>%
  mutate(pos = sprintf("%s_%s", chr, pos))

overlap.tmd450.tcga <- intersect(ann450k.hg38$pos,  cpg450df$cpg)
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(overlapTCGA = ifelse(cpg %in% overlap.tmd450.tcga, "yes", "no")) %>%
  mutate(name = ifelse(overlapTCGA == "yes", subset(ann450k.hg38, ann450k.hg38$pos == cpg)$Name, "none"))

if (file.exists(file.path(path.to.save.panel.info, "TMD450_overlapping_TCGA.xlsx")) == FALSE){
  writexl::write_xlsx(cpg450df, file.path(path.to.save.panel.info, "TMD450_overlapping_TCGA.xlsx"))  
}
  
##### Number of overlapped CpG in each region in TMD450
checkdf <- table(cpg450df$region, cpg450df$overlapTCGA) %>% as.data.frame() %>%
  subset(Var2 == "yes" & Freq != 0) %>% arrange(desc(Freq))

if (file.exists(file.path(path.to.save.panel.info, "count_CpG_in_regions_TMD450_overlapping_TCGA.xlsx")) == FALSE){
  writexl::write_xlsx(checkdf, file.path(path.to.save.panel.info, "count_CpG_in_regions_TMD450_overlapping_TCGA.xlsx"))  
}
