gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
source(file.path(path.to.main.src, "prepare_promoter_regions.R"))
library(GenomicRanges)
source(file.path(path.to.main.src, "configs.R"))

library(ggpubr)
# min.cov <- params$min.cov
# analysis.version <- params$analysis.version
# data.version <- params$data.version
# output.version <- params$output.version
# input.cancer.class <- params$input.cancer.class

# path.to.offline.pkgs <- "/media/hieunguyen/HNSD01/storage/offline_pkgs"
# install.packages(file.path(path.to.offline.pkgs, "methylKit_1.30.0.tar.gz"), type = "source", repos = NULL)
# install.packages(file.path(path.to.offline.pkgs, "BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz"), type = "source", repos = NULL)

analysis.version <- "0.1"
data.version <- "TMD_cov"
output.version <- "20240907"

min.cov.bases <- configs[[analysis.version]]$min.cov.bases
qvalue.cutoff <- configs[[analysis.version]]$qvalue.cutoff
methdiff.cutoff <- configs[[analysis.version]]$methdiff.cutoff
log2FC.cutoff <- configs[[analysis.version]]$log2FC.cutoff
up.flank <- configs[[analysis.version]]$up.flank
down.flank <- configs[[analysis.version]]$down.flank

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version)
path.to.save.QC.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version, "QC")

qcdf <- read.csv(file.path(path.to.save.QC.output, sprintf("metadata_cfDNA_lowpdepth_TMD_bam_cov.addCount450.%s.xlsx", 5)))
qcdf <- subset(qcdf, select = -c(num.CpG.tmd450))
for (i in c(5, 10, 15, 20)){
  file <- file.path(path.to.save.QC.output, sprintf("metadata_cfDNA_lowpdepth_TMD_bam_cov.addCount450.%s.xlsx", i))
  min.cov <- str_replace(basename(file), "metadata_cfDNA_lowpdepth_TMD_bam_cov.addCount450.", "")
  min.cov <- str_replace(min.cov, ".xlsx", "")
  tmpdf <- read.csv(file)[c("cov.name", "num.CpG.tmd450")]
  colnames(tmpdf) <- c("cov.name", sprintf("count.depth.%s", min.cov))
  qcdf <- merge(qcdf, tmpdf, by.x = "cov.name", by.y = "cov.name")
}

cpg450df <- read.csv(file.path(path.to.main.src, "all_cpg_in_450_tmd_regions.csv"))
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(cpg = sprintf("%s_%s", chrom, pos)) %>%
  subset(select = -c(X))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx"))

for (i in c(5, 10, 15, 20)){
  qcdf[[sprintf("pct.depth.%s", i)]] <- unlist(lapply(
    seq(1, nrow(qcdf)), function(x){
      qcdf[x, sprintf("count.depth.%s", i)]/qcdf[x, "num.CpG"]
    }
  ))
}

#####----------------------------------------------------------------------#####
##### total number of CpG in each samples
#####----------------------------------------------------------------------#####
my_comparisons <- list( c("Breast", "Control"), 
                        c("Liver", "Control"),
                        c("Lung", "Control"),
                        c("Gastric", "Control"),
                        c("CRC", "Control"))

qcdf %>% subset(select = c(num.CpG, cov.name, Label)) %>%
  ggplot(aes(x = Label,  y = num.CpG)) + geom_boxplot() + 
  stat_compare_means(comparisons = my_comparisons) + theme_pubr(base_size = 20) + theme(axis.text.x = element_text(angle = 90)) 

#####----------------------------------------------------------------------#####
##### total number of CpG in 450 regions at different depths
#####----------------------------------------------------------------------#####
for (i in c(5, 10, 15, 20)){
  p <- qcdf %>%
    ggplot(aes_string(x = "Label",  y = sprintf("count.depth.%s", i))) + geom_boxplot() + 
    stat_compare_means(comparisons = my_comparisons) + theme_pubr(base_size = 20) + theme(axis.text.x = element_text(angle = 90)) +
    ggtitle(sprintf("Number of CpG in 450 regions depth >= %s reads", i))  
}

#####----------------------------------------------------------------------#####
##### Number of CpG at each depth of coverage value
#####----------------------------------------------------------------------#####
plotdf <- qcdf[, c( "Label", to_vec( for(item in colnames(qcdf)) if (grepl("count.", item) == TRUE) item))] %>%
  pivot_longer(!Label, names_to = "group", values_to = "count")

plotdf$group <- factor(plotdf$group, levels = c("count.depth.5", 
                                                "count.depth.10",
                                                "count.depth.15",
                                                "count.depth.20"))

plotdf %>% ggplot(aes(x = group, y = count)) + geom_boxplot() + theme_pubr(base_size = 20) 

#####----------------------------------------------------------------------#####
##### on target rates
#####----------------------------------------------------------------------#####
plotdf <- qcdf %>%
  subset(Set == "train") 
plotdf <- plotdf[, c( "Label", to_vec( for(item in colnames(qcdf)) if (grepl("pct.", item) == TRUE) item))]  %>%
  pivot_longer(!Label, names_to = "group", values_to = "count")

plotdf$group <- factor(plotdf$group, levels = c("pct.depth.5", 
                                                "pct.depth.10",
                                                "pct.depth.15",
                                                "pct.depth.20"))

plotdf %>% ggplot(aes(x = group, y = count)) + geom_boxplot() + theme_pubr(base_size = 20) 

#####----------------------------------------------------------------------#####
##### check available CpG
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.save.QC.output, "count_occurrences_CpG_in_all_samples.rds")) == FALSE){
  all.cov.files <- Sys.glob(file.path(path.to.input, sprintf("filtered_%sreads_TMD450regions_cov", 5), "*.cov"))
  cpg450df <- cpg450df[!duplicated(cpg450df$cpg),]
  all.tmd.cpgs <- to_vec( for (item in seq(1, nrow(cpg450df))) 0)
  names(all.tmd.cpgs) <- cpg450df$cpg
  
  tmp.all.cov.files <- all.cov.files
  for (file in all.cov.files){
    print(length(setdiff(tmp.all.cov.files, c(file))))
    tmp.all.cov.files <- setdiff(tmp.all.cov.files, c(file))
    tmpdf <- read.csv(file, sep = "\t", header = FALSE)
    tmp.cpg <- tmpdf$V7
    check.count <- to_vec(for (i in names(all.tmd.cpgs)) if (i %in% tmp.cpg) 1 else 0)
    all.tmd.cpgs <- all.tmd.cpgs + check.count
  }
  saveRDS(all.tmd.cpgs, file.path(path.to.save.QC.output, "count_occurrences_CpG_in_all_samples.rds"))
  
} else {
  all.tmd.cpgs <- readRDS(file.path(path.to.save.QC.output, "count_occurrences_CpG_in_all_samples.rds"))
}


