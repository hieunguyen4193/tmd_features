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

##### panel overlapping CpG TMD450 and TCGA
path.to.save.panel.info <- file.path(path.to.main.output, "panel")
cpg450df <- readxl::read_excel(file.path(path.to.save.panel.info, "TMD450_overlapping_TCGA.xlsx")) %>%
  subset(overlapTCGA == "yes")
cpg450df <- cpg450df[!duplicated(cpg450df$cpg), ]

if (file.exists(file.path(path.to.04.output, sprintf("covdf_all_samples.rds"))) == FALSE){
  maindf <- data.frame(pos = cpg450df$cpg)
  tmp.all.cov.files <- all.cov.files
  for (file in all.cov.files){
    print(length(setdiff(tmp.all.cov.files, c(file))))
    tmp.all.cov.files <- setdiff(tmp.all.cov.files, c(file))
    tmpdf <- read.csv(file, sep = "\t", header = FALSE)[c("V7", "V4")]
    colnames(tmpdf) <- c("pos", str_replace(basename(file), ".cov", ""))
    maindf <- merge(maindf, tmpdf, by.x = "pos", by.y = "pos", all.x = TRUE)
  }
  saveRDS(maindf, file.path(path.to.04.output, sprintf("covdf_all_samples.rds")))
} else {
  maindf <- readRDS(file.path(path.to.04.output, sprintf("covdf_all_samples.rds")))
}

# colnames(maindf) <- to_vec( for(item in colnames(maindf)) str_replace_all(str_replace(item, "X", ""), "[.]", "-") )
maindf <- maindf %>% column_to_rownames("pos")

group1 <- subset(meta.data, meta.data$Label == input.cancer.class)$cov.name
group2 <- subset(meta.data, meta.data$Label == "Control")$cov.name

input.metadata <- data.frame(sample = c(group1, group2), 
                             label = c(
                               to_vec(for(item in seq(1, length(group1))) "cancer"),
                               to_vec(for(item in seq(1, length(group2))) "control")
                             ))

# this is the factor of interest
g <- factor(input.metadata$label, levels = c("cancer", "control"))
# use the above to create a design matrix
design <- model.matrix(~0+label, data=input.metadata)
colnames(design) <- levels(g)
fit <- lmFit(maindf[, intersect(colnames(maindf), input.metadata$sample)], design)
# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(cancer-control,
                            levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
resdf.raw <- topTable(fit2, num=Inf, coef=1)



