#' Read bismark coverage file as a methylKit object
#' 
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads coverage files,
#' which have chr,start,end, number of cytosines (methylated bases) 
#' and number of thymines (unmethylated bases).
#' 
#' @param location a list or vector of file paths to coverage files
#'     
#' @param sample.id a list or vector of sample ids
#' @param assembly a string for genome assembly. Any string would work.
#' @param treatment if there are multiple files to be read a treatment 
#'                  vector should be supplied.
#' @param context a string for context of methylation such as: "CpG" or "CHG"
#' @param min.cov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#' 
#' @return methylRaw or methylRawList objects
readBismarkCoverage<-function( location,sample.id,assembly="unknown",treatment,
                               context="CpG",min.cov=10)
{
  if(length(location)>1){
    stopifnot(length(location)==length(sample.id),
              length(location)==length(treatment))
  }
  
  result=list()
  for(i in 1:length(location)){
    df=fread.gzipped(location[[i]],data.table=FALSE)
    
    # remove low coverage stuff
    df=df[ (df[,5]+df[,6]) >= min.cov ,]
    
    
    
    
    # make the object (arrange columns of df), put it in a list
    result[[i]]= new("methylRaw",data.frame(chr=df[,1],start=df[,2],end=df[,3],
                                            strand="*",coverage=(df[,5]+df[,6]),
                                            numCs=df[,5],numTs=df[,6]),
                     sample.id=sample.id[[i]],
                     assembly=assembly,context=context,resolution="base"
    )
  }
  
  if(length(result) == 1){
    return(result[[1]])
  }else{
    
    new("methylRawList",result,treatment=treatment)
  }
  
}

#' Read bismark cytosine report file as a methylKit object
#' 
#' Bismark aligner can output methylation information per base in
#' multiple different formats. This function reads cytosine report files,
#' which have chr,start, strand, number of cytosines (methylated bases) 
#' and number of thymines (unmethylated bases),context, trinucletide context.
#' 
#' @param location a list or vector of file paths to coverage files
#'     
#' @param sample.id a list or vector of sample ids
#' @param assembly a string for genome assembly. Any string would work.
#' @param treatment if there are multiple files to be read a treatment 
#'                  vector should be supplied.
#' @param context a string for context of methylation such as: "CpG" or "CHG"
#' @param min.cov a numeric value for minimum coverage. Bases that have coverage
#' below this value will be removed.
#' 
#' @return methylRaw or methylRawList objects

readBismarkCytosineReport<-function(location,sample.id,assembly="unknown",treatment,
                                    context="CpG",min.cov=10){
  if(length(location)>1){
    stopifnot(length(location)==length(sample.id),
              length(location)==length(treatment))
  }
  
  result=list()
  for(i in 1:length(location)){
    df=fread.gzipped(location[[i]],data.table=FALSE)
    
    # remove low coverage stuff
    df=df[ (df[,4]+df[,5]) >= min.cov ,]
    
    
    
    
    # make the object (arrange columns of df), put it in a list
    result[[i]]= new("methylRaw",
                     data.frame(chr=df[,1],start=df[,2],end=df[,2],
                                strand=df[,3],coverage=(df[,4]+df[,5]),
                                numCs=df[,4],numTs=df[,5]),
                     sample.id=sample.id[[i]],
                     assembly=assembly,context=context,resolution="base"
    )
  }
  
  if(length(result) == 1){
    return(result[[1]])
  }else{
    
    new("methylRawList",result,treatment=treatment)
  }
}

# reads gzipped files,
fread.gzipped<-function(filepath,...){
  require(R.utils)
  require(data.table)
  
  
  
  # decompress first, fread can't read gzipped files
  if (R.utils::isGzipped(filepath)){
    
    if(.Platform$OS.type == "unix") {
      filepath=paste("zcat",filepath)
    } else {
      filepath <- R.utils::gunzip(filepath,temporary = FALSE, overwrite = TRUE,
                                  remove = FALSE)
    }
    
    
  }
  
  ## Read in the file
  fread(filepath,...)
  
}

`%ni%` = Negate(`%in%`)

# Some SMALL helper functions
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

annotate_with_promoter_database <- function(input.grange, ref.grange, add_suffix = "promoter"){
  overlap.idxs <-  findOverlaps(input.grange, ref.grange)
  tmpdf1 <- input.grange[queryHits(overlap.idxs), ] %>% data.frame() %>%
    rowwise() %>%
    mutate(name = sprintf("%s.%s.%s", seqnames, start, end))
  tmpdf2 <- ref.grange[subjectHits(overlap.idxs), ] %>% data.frame() %>%
    rowwise() %>%
    mutate(promoter.name = sprintf("%s.%s.%s", seqnames, start, end))
  colnames(tmpdf2) <- to_vec ( for (item in colnames(tmpdf2)) sprintf("%s_%s", add_suffix, item))
  outputdf <- cbind(subset(tmpdf1, select = c(name)),
                    tmpdf2)
  return(outputdf)
}
