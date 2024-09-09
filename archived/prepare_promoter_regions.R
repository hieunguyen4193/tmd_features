#> A database containing all human-genome TSS was downloaded from UCSC.
#> Based on this database, we define a region of x bp upstream and y bp downstream
#> as promoter. We intersect our DMR/DML with these promoters.
#> Refseq database was also downloaded from UCSC, choose refseq ALL.
#> -----

define_promoter_regions <- function(up.flank, down.flank, path.to.save.promoterdf){
  if (file.exists(file.path(path.to.save.promoterdf, sprintf("promoter_regions_up_%s_down_%s_from_UCSC_TSS.csv", up.flank, down.flank))) == FALSE){
    library(GenomicRanges)
    library(dplyr)
    library(tidyverse)
    library(comprehenr)
    
    path.to.main.src.tmp <- "/media/hieunguyen/HNSD01/src/PBMC"
    path.to.TSS.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.bed")
    path.to.refseq <- file.path(path.to.main.src.tmp, "hg19_refseq_all.bed")
    
    path.to.TSS.annotated.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.annotated.bed")
    
    hg19.refseq <- read.csv(file.path(path.to.refseq), sep = "\t")[,c("chrom", "txStart", "txEnd", "name", "name2", "strand")]
    all.valid.chroms <- to_vec( for(item in seq(1,22)) sprintf("chr%s", item))
    hg19.refseq <- subset(hg19.refseq, hg19.refseq$chrom %in% all.valid.chroms)
    colnames(hg19.refseq) <- c("chrom", "start", "end", "tx_name", "gene", "strand")
    hg19.refseq.grange <- makeGRangesFromDataFrame(df = hg19.refseq, start.field = "start", end.field = "end", seqnames.field = "chrom", keep.extra.columns = TRUE, strand.field = "strand")
    
    tssdf <- read.csv(path.to.TSS.file, sep = "\t", header = FALSE)[c("V2", "V3", "V4", "V7", "V12")]
    colnames(tssdf) <- c("chrom", "start", "end", "strand", "pseudogene")
    
    tss.grange <- makeGRangesFromDataFrame(df = tssdf, 
                                           start.field = "start", 
                                           end.field = 'end', 
                                           seqnames.field = "chrom", 
                                           keep.extra.columns = TRUE)
    
    
    if (file.exists(path.to.TSS.annotated.file) == FALSE){
      add_suffix1 <- "TSS"
      add_suffix2 <- "hg19"
      
      overlap.idxs <-  findOverlaps(tss.grange, hg19.refseq.grange)
      tmpdf1 <- tss.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
        rowwise() %>%
        mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
      colnames(tmpdf1) <- to_vec ( for (item in colnames(tmpdf1)) sprintf("%s_%s", add_suffix1, item))
      tmpdf2 <- hg19.refseq.grange[subjectHits(overlap.idxs), ] %>% data.frame() %>%
        rowwise() %>%
        mutate(promoter.name = sprintf("%s.%s.%s", seqnames, start, end))
      colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix2, item))
      
      tss.annotated.df <- cbind(tmpdf1, tmpdf2)
      write.csv(tss.annotated.df, path.to.TSS.annotated.file)
    }
    
    promoterdf <- tssdf %>% rowwise() %>% 
      mutate(promoter.start = ifelse(strand == "+", start - up.flank, end - down.flank)) %>%
      mutate(promoter.end = ifelse(strand == "+", start + down.flank, end + up.flank))
    
    promoter.grange <- makeGRangesFromDataFrame(df = subset(promoterdf, select = c(chrom, promoter.start, promoter.end, strand)), 
                                                seqnames.field = "chrom",
                                                start.field = "promoter.start",
                                                end.field = "promoter.end",
                                                strand.field = "strand",
                                                keep.extra.columns = TRUE)
    
    add_suffix1 <- "promoter"
    add_suffix2 <- "hg19"
    
    overlap.idxs <-  findOverlaps(promoter.grange, hg19.refseq.grange)
    
    tmpdf1 <- promoter.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
      rowwise() %>%
      mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
    colnames(tmpdf1) <- to_vec ( for (item in colnames(tmpdf1)) sprintf("%s_%s", add_suffix1, item))
    
    tmpdf2 <- hg19.refseq.grange[subjectHits(overlap.idxs), ] %>% data.frame() 
    colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix2, item))
    
    promoterdf <- cbind(tmpdf1, tmpdf2)
    write.csv(promoterdf, file.path(path.to.save.promoterdf, sprintf("promoter_regions_up_%s_down_%s_from_UCSC_TSS.csv", up.flank, down.flank)))
  } else {
    promoterdf <- read.csv(file.path(path.to.save.promoterdf, sprintf("promoter_regions_up_%s_down_%s_from_UCSC_TSS.csv", up.flank, down.flank))) %>%
      subset(select = -c(X))
  }
  
  return(promoterdf)
}

##### TEST 
library(GenomicRanges)
library(dplyr)
library(tidyverse)
library(comprehenr)

path.to.main.src.tmp <- "/media/hieunguyen/HNSD01/src/PBMC"
path.to.TSS.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.bed")
path.to.refseq <- file.path(path.to.main.src.tmp, "hg19_refseq_all.bed")

path.to.TSS.annotated.file <- file.path(path.to.main.src.tmp, "UCSC_TSS.annotated.bed")
tss.annotated.df <- read.csv(path.to.TSS.annotated.file)
