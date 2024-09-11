gc()
rm(list = ls())

#> DESCRIPTION
#> We count the number of Cpg sites in the TMD 450 regions (on target) and 
#> off target CpG sites. See how many CpG a sample has. 
#> 
path.to.main.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))

data.version <- "TMD_cov"
output.version <- "20240910"
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "TMD450_TCGA_data_analysis"
path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, PROJECT, output.version)
path.to.save.panel.info <- file.path(path.to.main.output, "panel")

for (min.cov in c(5, 10, 15, 20)){
  all.cov.files <- Sys.glob(file.path(path.to.input, sprintf("filtered_%sreads_cov", min.cov), "*.cov"))
  
  path.to.save.filtered450.cov <- file.path(path.to.input, sprintf("filtered_%sreads_TMD450regionsOverlapTCGA_cov", min.cov))
  dir.create(path.to.save.filtered450.cov, showWarnings = FALSE, recursive = TRUE)
  
  names(all.cov.files) <- unlist(lapply(all.cov.files, function(x){
    x <- basename(x)
    x <- str_split(x, ".deduplicated")[[1]][[1]]
    return(x)
  }))
  
  cpg450df <- readxl::read_excel(file.path(path.to.save.panel.info, "TMD450_overlapping_TCGA.xlsx"))
  cpg450df <- subset(cpg450df, cpg450df$overlapTCGA == "yes")
  
  for (sample.id in names(all.cov.files)){
    print(sprintf("working on sample %s", sample.id))
    if (file.exists(file.path(path.to.save.filtered450.cov, sprintf("%s.cov", sample.id))) == FALSE){
      tmpdf <- read.csv(all.cov.files[[sample.id]], sep = "\t", header = FALSE)
      tmpdf <- tmpdf %>% rowwise() %>%
        mutate(cpg = sprintf("%s_%s", V1, V2))
      tmpdf <- subset(tmpdf, tmpdf$cpg %in% cpg450df$cpg)
      write.table(tmpdf, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE,  
                  file = file.path(path.to.save.filtered450.cov, sprintf("%s.cov", sample.id)))
    } else {
      if (file.size(file.path(path.to.save.filtered450.cov, sprintf("%s.cov", sample.id))) == 0L){
        tmpdf <- data.frame()
      } else {
        tmpdf <- read.csv(file.path(path.to.save.filtered450.cov, sprintf("%s.cov", sample.id)), sep = "\t", header = FALSE, )      
      }
    }
  }
  
}

