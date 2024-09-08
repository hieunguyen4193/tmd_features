# gc()
# rm(list = ls())

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

min.cov <- 5
analysis.version <- "0.1"
data.version <- "TMD_cov"
output.version <- "20240907_sampling_Control"
input.cancer.class <- "CRC"
all.cancer.classes <- c("Liver", "Breast", "Gastric", "Lung", "CRC")

min.cov.bases <- configs[[analysis.version]]$min.cov.bases
qvalue.cutoff <- configs[[analysis.version]]$qvalue.cutoff
methdiff.cutoff <- configs[[analysis.version]]$methdiff.cutoff
log2FC.cutoff <- configs[[analysis.version]]$log2FC.cutoff
up.flank <- configs[[analysis.version]]$up.flank
down.flank <- configs[[analysis.version]]$down.flank

outdir <- "/media/hieunguyen/HNSD_mini/outdir"
path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, "TMD_read_based_features", "output", sprintf("data_%s", data.version), output.version)

path.to.01.output <- file.path(path.to.main.output, sprintf("01_output_%s", analysis.version), sprintf("minCov_%s_class_%s", min.cov, input.cancer.class))
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

all.cov.files <- Sys.glob(file.path(path.to.input, sprintf("filtered_%sreads_cov", min.cov), "*.cov"))
names(all.cov.files) <- unlist(lapply(all.cov.files, function(x){
  x <- basename(x)
  x <- str_split(x, ".deduplicated")[[1]][[1]]
  return(x)
}))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "metadata_cfDNA_lowpdepth_TMD_bam_cov.xlsx"))
meta.data <- subset(meta.data, meta.data$Label %in% c("Control", input.cancer.class))
meta.data <- meta.data %>% rowwise() %>%
  mutate(cov.name = str_split(basename(TM_COV), ".deduplicated")[[1]][[1]])
meta.data <- subset(meta.data, is.na(meta.data$cov.name) == FALSE)
meta.data <- meta.data[!duplicated(meta.data$cov.name), ]

count.CpG <- read.csv("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov/TMD_cov/count_CpG_filtered_5reads_cov.txt", sep = "\t", header = FALSE)
count.CpG <- count.CpG %>% rowwise() %>%
  mutate(cov.name = str_split(basename(V1), ".deduplicated")[[1]][[1]])

meta.data.control <- subset(meta.data, meta.data$Label == "Control")
meta.data.cancer <- subset(meta.data, meta.data$Label == input.cancer.class)

n_cancer <- nrow(meta.data.cancer)
n_control <- nrow(meta.data.control)

if (file.exists(file.path(path.to.01.output, "sampling.control.samples.rds")) == FALSE){
  set.seed(411)
  sampling.control.samples <- sample(meta.data.control$cov.name, nrow(meta.data.cancer))
  saveRDS(sampling.control.samples, file.path(path.to.01.output, "sampling.control.samples.rds"))  
} else {
  sampling.control.samples <- readRDS(file.path(path.to.01.output, "sampling.control.samples.rds"))
}

meta.data <- subset(meta.data, meta.data$cov.name %in% c(meta.data.cancer$cov.name, sampling.control.samples))

all.cov.files <- all.cov.files[meta.data$cov.name]
labels <- to_vec( for(item in names(all.cov.files)) if (subset(meta.data, meta.data$cov.name == item)$Label == input.cancer.class) 1 else 0)

#####----------------------------------------------------------------------#####
# Generate DML object
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "DML_obj.rds")) == FALSE){
  print("Generate methylkit object from input cov files...")
  DML.obj <- readBismarkCoverage( all.cov.files,
                                  sample.id = names(all.cov.files),
                                  assembly = "hg38",
                                  treatment = labels,
                                  context = "CpG",
                                  min.cov = min.cov)
  saveRDS(DML.obj, file.path(path.to.01.output, "DML_obj.rds"))  
} else {
  print("Methylkit object exists, reading in...")
  DML.obj <-  readRDS(file.path(path.to.01.output, "DML_obj.rds"))
}


#####----------------------------------------------------------------------#####
# Unite the meth objects
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "meth.rds")) == FALSE){
  print("Unite the methylation methylkit object...")
  meth <- methylKit::unite(object = DML.obj, 
                           destrand = FALSE, mc.cores = 40)
  saveRDS(meth, file.path(path.to.01.output, "meth.rds"))  
} else {
  print("Methylkit object united, reading in...")
  meth <- readRDS(file.path(path.to.01.output, "meth.rds"))
}

num.cpgs <- c()
all.samples <- c()
for (i in seq(1,length(DML.obj))){
  num.cpgs <- c(num.cpgs, nrow(getData(DML.obj[[i]])))
  all.samples <- c(all.samples, DML.obj[[i]]@sample.id)
}

countdf <- data.frame(cov.name = all.samples, num = num.cpgs)
countdf <- merge(countdf, meta.data, by.x = "cov.name", by.y = "cov.name")
countdf <- countdf %>% arrange(desc(num))
countdf %>% ggplot(aes(x = Label, y = num)) + geom_boxplot()
#####----------------------------------------------------------------------#####
# Generate PCA from CpG loci
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "pca_from_CpG_loci_res.rds")) == FALSE){
  print("Generate PCA from CpG methylation loci information...")
  pca.res <- PCASamples(meth, obj.return = TRUE, screeplot = FALSE)
} else {
  print("PCA done, reading in...")
  pca.res <- readRDS(file.path(path.to.01.output, "pca_from_CpG_loci_res.rds"))
}
pcadf <- pca.res$x %>% as.data.frame() %>% rownames_to_column("SampleID")  
pcadf <- merge(pcadf, meta.data, by.x = "SampleID", by.y = "cov.name")
pcadf %>% ggplot(aes(x = PC1, y = PC2, color = Label)) + geom_point(size = 3) + scale_color_manual(values = c("gray", "red"))

#####----------------------------------------------------------------------#####
# Generate DML
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.01.output, "DML.rds")) == FALSE){
  print("Conduct differential methylated test between case and control ...")
  myDiff <- calculateDiffMeth(meth, mc.cores = 45)
  saveRDS(myDiff, file.path(path.to.01.output, "DML.rds"))
} else {
  print("Differential methylated test results existed, reading in...")
  myDiff <- readRDS(file.path(path.to.01.output, "DML.rds"))
}

# and get all diff loci
if (file.exists(file.path(path.to.01.output, "diff_locidf.rds")) == FALSE){
  diff.loci <- getMethylDiff(myDiff, difference = methdiff.cutoff, qvalue = qvalue.cutoff)
  saveRDS(diff.loci, file.path(path.to.01.output, "diff_locidf.rds"))  
} else {
  diff.loci <- readRDS(file.path(path.to.01.output, "diff_locidf.rds"))
}

diff.locidf <- getData(diff.loci)

