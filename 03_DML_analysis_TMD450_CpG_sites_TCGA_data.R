gc()
rm(list = ls())

#> DESCRIPTION
#> We count the number of Cpg sites in the TMD 450 regions (on target) and 
#> off target CpG sites. See how many CpG a sample has. 
#> 
path.to.main.src <- "/media/hieunguyen/HNSD01/src/tmd_features"
source(file.path(path.to.main.src, "import_libraries.R"))
source(file.path(path.to.main.src, "helper_functions.R"))
library(limma)

##### input args
data.version <- "TMD_cov"
output.version <- "20240910"
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "TMD450_TCGA_data_analysis"

all.cancer.classes <- c("Liver", "Lung", "Breast", "CRC", "Gastric")

path.to.input <- file.path("/media/hieunguyen/GSHD_HN01/raw_data/bismark_cov", data.version)
path.to.main.output <- file.path(outdir, PROJECT, output.version)

for (input.cancer.class in all.cancer.classes){
# for (input.cancer.class in c("Liver")){
  path.to.01.output <- file.path(path.to.main.output, "01_output", input.cancer.class)
  path.to.03.output <- file.path(path.to.main.output, "03_output", input.cancer.class)
  dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)
  
  path.to.save.panel.info <- file.path(path.to.main.output, "panel")
  if (file.exists(file.path(path.to.03.output, "DMPs.xlsx")) == FALSE){
    ##### read TCGA idat preprocessed object
    maindf <- readRDS(file.path(path.to.01.output, "bVals.rds"))
    
    ##### panel overlapping CpG TMD450 and TCGA
    cpg450df <- readxl::read_excel(file.path(path.to.save.panel.info, "TMD450_overlapping_TCGA.xlsx")) %>%
      subset(overlapTCGA == "yes")
    cpg450df <- cpg450df[!duplicated(cpg450df$cpg), ]
    maindf <- maindf[intersect(row.names(maindf), cpg450df$name), ]
    
    ##### metadata
    path.to.tcga.info <- file.path(path.to.main.src, "TCGA_database_info")
    if (input.cancer.class == "CRC"){
      sample.sheet.normal <- read.csv(Sys.glob(file.path(path.to.tcga.info, sprintf("%s_idat", "Colon"), "*Normal_450K.tsv")), sep = "\t")
      sample.sheet.cancer <- read.csv(Sys.glob(file.path(path.to.tcga.info, sprintf("%s_idat", "Colon"), "*Tumor_450K.tsv")), sep = "\t")
    } else {
      sample.sheet.normal <- read.csv(Sys.glob(file.path(path.to.tcga.info, sprintf("%s_idat", input.cancer.class), "*Normal_450K.tsv")), sep = "\t")
      sample.sheet.cancer <- read.csv(Sys.glob(file.path(path.to.tcga.info, sprintf("%s_idat", input.cancer.class), "*Tumor_450K.tsv")), sep = "\t")
    }
    
    sample.sheet <- rbind(sample.sheet.normal, sample.sheet.cancer)
    sample.sheet <- sample.sheet %>% rowwise() %>%
      mutate(SampleID = str_replace(str_replace(File.Name, "_Grn.idat", ""), "_Red.idat", "")) %>%
      mutate(Label = ifelse(Sample.Type == "Solid Tissue Normal", "Control", "Tumor"))
    sample.sheet <- sample.sheet[!duplicated(sample.sheet$SampleID),]
    sample.sheet <- subset(sample.sheet, select = c(SampleID, Label))
    
    ##### run DML
    group1 <- subset(sample.sheet, sample.sheet$Label == "Tumor")$SampleID
    group2 <- subset(sample.sheet, sample.sheet$Label == "Control")$SampleID
    
    input.metadata <- data.frame(sample = c(group1, group2), 
                                 label = c(
                                   to_vec(for(item in seq(1, length(group1))) "Tumor"),
                                   to_vec(for(item in seq(1, length(group2))) "Control")
                                 ))
    g <- factor(input.metadata$label, levels = c("Tumor", "Control"))
    design <- model.matrix(~0+label, data=input.metadata)
    colnames(design) <- levels(g)
    fit <- lmFit(maindf[, input.metadata$sample], design)
    contMatrix <- makeContrasts(Tumor-Control,
                                levels=design)
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)
    
    # >>>> diff = cancer - control --> hyper cancer > control, hypo cancer < control
    diff.methyl <- rowMeans(maindf[, group1]) - rowMeans(maindf[, group2]) %>% as.data.frame() 
    
    diff.methyl <- diff.methyl %>% rownames_to_column("cpg") 
    colnames(diff.methyl) <- c("cpg", "diff")
    
    DMPs <- topTable(fit2, num=Inf, coef=1) %>% as.data.frame() %>% arrange(desc(logFC)) %>%
      rownames_to_column("TCGA_CpG") %>%
      rowwise() %>%
      mutate(cpg = subset(cpg450df, cpg450df$name == TCGA_CpG)$cpg) %>%
      mutate(hyper_or_hypo1 = ifelse(logFC >= 0, "hypo", "hyper"))
    
    DMPs <- merge(DMPs, diff.methyl, by.x = "TCGA_CpG", by.y = "cpg") %>%
      subset(adj.P.Val <= 0.05) %>%
      rowwise() %>%
      mutate(abs.diff = abs(diff)) 
    DMPs <- merge(DMPs, subset(cpg450df, select = c(cpg, region)), by.x = "cpg", by.y = "cpg")
    writexl::write_xlsx(DMPs, file.path(path.to.03.output, "DMPs.xlsx"))
    
    countDMPs <- table(DMPs$region, DMPs$hyper_or_hypo1) %>% as.data.frame() %>%
      pivot_wider(names_from = "Var2", values_from = "Freq")
    writexl::write_xlsx(countDMPs, file.path(path.to.03.output, "countDMPs.xlsx"))
  }
}
