gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
source(file.path(path.to.main.src, "prepare_promoter_regions.R"))
library(GenomicRanges)
source(file.path(path.to.main.src, "configs.R"))
library(limma)
library(ggpubr)

data.version <- "TMD_cov"
output.version <- "20240907"

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version)
path.to.difftest.output <- file.path(path.to.main.output, "difftest_output")
dir.create(path.to.difftest.output, showWarnings = FALSE, recursive = TRUE)

path.to.save.QC.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version, "QC")

all.cov.files <- Sys.glob(file.path(path.to.input, sprintf("filtered_%sreads_TMD450regions_cov", 5), "*.cov"))

##### get the summary counts of all 450 TMD regions CpGs
all.tmd.cpgs <- readRDS(file.path(path.to.save.QC.output, "count_occurrences_CpG_in_all_samples.rds"))
thres.CpG.occurrence <- round(median(all.tmd.cpgs))

##### read main df
if (file.exists(file.path(path.to.difftest.output, sprintf("covdf_all_samples.%s.csv", thres.CpG.occurrence))) == FALSE){
  maindf <- data.frame(pos = selected.CpGs)
  tmp.all.cov.files <- all.cov.files
  for (file in all.cov.files){
    print(length(setdiff(tmp.all.cov.files, c(file))))
    tmp.all.cov.files <- setdiff(tmp.all.cov.files, c(file))
    tmpdf <- read.csv(file, sep = "\t", header = FALSE)[c("V7", "V4")]
    colnames(tmpdf) <- c("pos", str_replace(basename(file), ".cov", ""))
    maindf <- merge(maindf, tmpdf, by.x = "pos", by.y = "pos", all.x = TRUE)
  }
  write.csv(maindf, file.path(path.to.difftest.output, sprintf("covdf_all_samples.%s.csv", thres.CpG.occurrence)))
} else {
  maindf <- read.csv(file.path(path.to.difftest.output, sprintf("covdf_all_samples.%s.csv", thres.CpG.occurrence)))
  maindf <- subset(maindf, select = -c(X))
}
colnames(maindf) <- to_vec( for(item in colnames(maindf)) str_replace_all(str_replace(item, "X", ""), "[.]", "-") )
maindf <- maindf %>% column_to_rownames("pos")

##### get the metadata
meta.data <- readxl::read_excel(file.path(path.to.main.src, "metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx"))

meta.data <- meta.data %>% rowwise() %>%
  mutate(cov.name = str_split(basename(TM_COV), ".deduplicated")[[1]][[1]])
meta.data <- subset(meta.data, is.na(meta.data$cov.name) == FALSE)
meta.data <- subset(meta.data, meta.data$cov.name %in% colnames(maindf))

meta.data.train <- subset(meta.data, meta.data$Set == "train")
  
train.samples <- list()
for (input.class in unique(meta.data.train$Label)){
  train.samples[[input.class]] <- subset(meta.data.train, meta.data.train$Label == input.class)$cov.name
}


selected.CpGs <- names(all.tmd.cpgs[all.tmd.cpgs >= thres.CpG.occurrence])

#####----------------------------------------------------------------------#####
##### get all CpG information from the TMD 450 regions
#####----------------------------------------------------------------------#####
cpg450df <- read.csv(file.path(path.to.main.src, "all_cpg_in_450_tmd_regions.csv"))
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(cpg = sprintf("%s_%s", chrom, pos)) %>%
  subset(select = -c(X))

#####----------------------------------------------------------------------#####
##### Run differential methylated test for each CpG
#####----------------------------------------------------------------------#####
all.check.region <- list()

for (group in c("Liver", "Breast", "CRC", "Lung", "Gastric")){
  group1 <- train.samples[[group]]
  group2 <- train.samples[["Control"]]
  
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
  fit <- lmFit(maindf[, input.metadata$sample], design)
  # create a contrast matrix for specific comparisons
  contMatrix <- makeContrasts(cancer-control,
                              levels=design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  resdf.raw <- topTable(fit2, num=Inf, coef=1)
  resdf <- subset(resdf.raw, resdf.raw$adj.P.Val <= 0.05) %>% rownames_to_column("cpg") %>%
    rowwise() %>%
    mutate(hypo_or_hyper = ifelse(logFC > 0, "hyper", "hypo"))
  
  resdf <- merge(resdf, subset(cpg450df, select = c(region, cpg)), by.x = "cpg", by.y = "cpg")
  
  check.region.longer <- table(resdf$region, resdf$hypo_or_hyper) %>% data.frame() 
  check.region <- check.region.longer %>% pivot_wider(names_from = "Var2", values_from = "Freq")
  colnames(check.region) <- c("CpG", "hyper", "hypo")
  
  check.region <- check.region %>% rowwise() %>%
    mutate(diff.count = hyper - hypo) %>%
    mutate(region_type = ifelse(hyper > hypo, "hyper", "hypo"))
  
  writexl::write_xlsx(resdf, file.path(path.to.difftest.output, sprintf("diff_test_%s.xlsx", group)))
  writexl::write_xlsx(check.region, file.path(path.to.difftest.output, sprintf("diff_test_%s_region_summary.xlsx", group)))
  all.check.region[[group]] <- check.region
  
  print(group)
  print(table(all.check.region[[group]]$region_type))
}

