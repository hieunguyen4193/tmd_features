gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
source(file.path(path.to.main.src, "prepare_promoter_regions.R"))
library(GenomicRanges)
source(file.path(path.to.main.src, "configs.R"))

# min.cov <- params$min.cov
# analysis.version <- params$analysis.version
# data.version <- params$data.version
# output.version <- params$output.version
# input.cancer.class <- params$input.cancer.class

# path.to.offline.pkgs <- "/media/hieunguyen/HNSD01/storage/offline_pkgs"
# install.packages(file.path(path.to.offline.pkgs, "methylKit_1.30.0.tar.gz"), type = "source", repos = NULL)
# install.packages(file.path(path.to.offline.pkgs, "BSgenome.Hsapiens.UCSC.hg38_1.4.3.tar.gz"), type = "source", repos = NULL)

min.cov <- 15
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
dir.create(path.to.save.QC.output, showWarnings = FALSE, recursive = TRUE)

all.cov.files <- Sys.glob(file.path(path.to.input, sprintf("filtered_%sreads_cov", min.cov), "*.cov"))

path.to.save.filtered450.cov <- file.path(path.to.input, sprintf("filtered_%sreads_TMD450regions_cov", min.cov))
dir.create(path.to.save.filtered450.cov, showWarnings = FALSE, recursive = TRUE)

names(all.cov.files) <- unlist(lapply(all.cov.files, function(x){
  x <- basename(x)
  x <- str_split(x, ".deduplicated")[[1]][[1]]
  return(x)
}))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx"))
meta.data <- meta.data %>% rowwise() %>%
  mutate(cov.name = str_split(basename(TM_COV), ".deduplicated")[[1]][[1]])
meta.data <- subset(meta.data, is.na(meta.data$cov.name) == FALSE)
meta.data <- meta.data[!duplicated(meta.data$cov.name), ]

count.CpG <- read.csv("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov/TMD_cov/count_CpG_filtered_5reads_cov.txt", sep = "\t", header = FALSE)
count.CpG <- count.CpG %>% rowwise() %>%
  mutate(cov.name = str_split(basename(V1), ".deduplicated")[[1]][[1]])
colnames(count.CpG) <- c("filename", "num.CpG", "cov.name")

meta.data <- merge(meta.data, subset(count.CpG, select = c(cov.name, num.CpG)), by.x = "cov.name", by.y = "cov.name")

meta.data %>% ggplot(aes(x = Label, y = num.CpG, fill = Label)) + 
  geom_boxplot() +
  theme(axis.text = element_text(angle = 90))

cpg450df <- read.csv(file.path(path.to.main.src, "all_cpg_in_450_tmd_regions.csv"))
cpg450df <- cpg450df %>% rowwise() %>%
  mutate(cpg = sprintf("%s_%s", chrom, pos))

all.samples <- c()
tmp.counts <- c()
for (sample.id in names(all.cov.files)){
  print(sprintf("working on sample %s", sample.id))
  if (file.exists(file.path(path.to.save.filtered450.cov, sprintf("%s.tsv", sample.id))) == FALSE){
    tmpdf <- read.csv(all.cov.files[[sample.id]], sep = "\t", header = FALSE)
    tmpdf <- tmpdf %>% rowwise() %>%
      mutate(cpg = sprintf("%s_%s", V1, V2))
    tmpdf <- subset(tmpdf, tmpdf$cpg %in% cpg450df$cpg)
    write.table(tmpdf, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE,  
                file = file.path(path.to.save.filtered450.cov, sprintf("%s.tsv", sample.id)))
  } else {
    if (file.size(file.path(path.to.save.filtered450.cov, sprintf("%s.tsv", sample.id))) == 0L){
      tmpdf <- data.frame()
    } else {
      tmpdf <- read.csv(file.path(path.to.save.filtered450.cov, sprintf("%s.tsv", sample.id)), sep = "\t", header = FALSE, )      
    }

  }
  tmp.counts <- c(tmp.counts, nrow(tmpdf))
  all.samples <- c(all.samples, sample.id)   
}

count.CpG.in.tmd450 <- data.frame(cov.name = all.samples, num.CpG.tmd450 = tmp.counts)
meta.data <- merge(meta.data, count.CpG.in.tmd450, by.x = "cov.name", by.y = "cov.name")
write.csv(meta.data, file.path(path.to.save.QC.output, sprintf("metadata_cfDNA_lowpdepth_TMD_bam_cov.addCount450.%s.xlsx", min.cov)))
