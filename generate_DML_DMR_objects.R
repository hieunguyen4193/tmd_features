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
# install.packages(file.path(path.to.offline.pkgs, "BSgenome.Hsapiens.UCSC.hg19_1.4.3.tar.gz"), type = "source", repos = NULL)

min.cov <- 5
analysis.version <- "0.1"
data.version <- "TMD_cov"
output.version <- "20240907"
# input.cancer.class <- "Liver"
all.cancer.classes <- c("Liver", "Breast", "Gastric", "Lung", "CRC")

for (input.cancer.class in all.cancer.classes){
  
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
  all.cov.files <- all.cov.files[meta.data$cov.name]
  labels <- to_vec( for(item in names(all.cov.files)) if (subset(meta.data, meta.data$cov.name == item)$Label == input.cancer.class) 1 else 0)
  #####----------------------------------------------------------------------#####
  # Generate DML object
  #####----------------------------------------------------------------------#####
  
  if (file.exists(file.path(path.to.01.output, "DML_obj.rds")) == FALSE){
    print("Generate methylkit object from input cov files...")
    DML.obj <- readBismarkCoverage( all.cov.files,
                                    sample.id = names(all.cov.files),
                                    assembly = "hg19",
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
    meth <- methylKit::unite(object = DML.obj, destrand = FALSE)
    saveRDS(meth, file.path(path.to.01.output, "meth.rds"))  
  } else {
    print("Methylkit object united, reading in...")
    meth <- readRDS(file.path(path.to.01.output, "meth.rds"))
  }
  
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
  
  #####----------------------------------------------------------------------#####
  # Prepare data for DMR analysis
  #####----------------------------------------------------------------------#####
  methdf <- getData(meth)
  methdf.grange <- makeGRangesFromDataFrame(methdf)
  
  if (file.exists(file.path(path.to.01.output, "finished_generating_DMR.csv")) == FALSE){
    print("Prepare data for DMR analysis...")
    meth.chr <- hash()
    tiles <- hash()
    diff.regions.by.chr <- hash()
    methdf.grange.chr <- hash()
    
    diff.regiondf <- data.frame()
    for (chrom in setdiff(unique(methdf$chr), c("X", "Y", "MT"))){
      print(sprintf("Working on chromosome: %s", chrom))
      methdf.grange.chr[[chrom]] <- subset(methdf.grange, seqnames == chrom) 
      
      meth.chr[[chrom]] <- selectByOverlap(meth, methdf.grange.chr[[chrom]])
      meth.chr[[chrom]] <- meth.chr[[chrom]][!duplicated(meth.chr[[chrom]]), ] 
      
      tiles[[chrom]] <- tileMethylCounts(meth.chr[[chrom]], win.size=1000, step.size=1000, cov.bases = min.cov.bases, mc.cores = 40)
      
      if( dim(tiles[[chrom]])[[1]] != 0){
        diff.regions.by.chr[[chrom]] <- data.frame(calculateDiffMeth(tiles[[chrom]])) 
        diff.regiondf <- rbind(diff.regiondf, diff.regions.by.chr[[chrom]])
      }
    }
    if (nrow(diff.regiondf) != 0){
      ##### diff meth cut-off not applied in diff regions
      diff.regiondf <- diff.regiondf %>% rowwise() %>%
        mutate(abs.meth.diff = abs(meth.diff)) %>%
        mutate(sig = ifelse(qvalue <= qvalue.cutoff, "significant", "not.significant"))
      diff.regiondf.raw <- diff.regiondf
      
      diff.regiondf <- subset(diff.regiondf, sig == "significant")
      diff.regiondf <- diff.regiondf %>% rowwise() %>% 
        mutate(name = sprintf("%s.%s.%s", chr,  start, end))
      
      saveRDS(meth.chr, file.path(path.to.01.output, "meth_chr.rds"))
      saveRDS(tiles, file.path(path.to.01.output, "tiles.rds"))
      saveRDS(diff.regions.by.chr, file.path(path.to.01.output, "diff_regions_by_chr.rds"))
      saveRDS(methdf.grange.chr, file.path(path.to.01.output, "methdf_grange_chr.rds"))
      saveRDS(diff.regiondf.raw, file.path(path.to.01.output, "diff_regiondf_raw.rds"))
      saveRDS(diff.regiondf, file.path(path.to.01.output, "diff_regiondf.rds"))
      write.csv(data.frame(status = c("Finished generating DMR results")), file.path(path.to.01.output, "finished_generating_DMR.csv"))    
    }
  } else {
    if (file.exists(file.path(path.to.01.output, "diff_regiondf.rds")) == TRUE){
      print("Reading in DMR information...")
      meth.chr <- readRDS(file.path(path.to.01.output, "meth_chr.rds"))
      tiles <- readRDS(file.path(path.to.01.output, "tiles.rds"))
      diff.regions.by.chr <- readRDS(file.path(path.to.01.output, "diff_regions_by_chr.rds"))
      methdf.grange.chr <- readRDS(file.path(path.to.01.output, "methdf_grange_chr.rds"))
      diff.regiondf.raw <- readRDS(file.path(path.to.01.output, "diff_regiondf_raw.rds"))
      diff.regiondf <- readRDS(file.path(path.to.01.output, "diff_regiondf.rds"))
    } else {
      print("cannot perform DMR analysis")
    }
  }
}

