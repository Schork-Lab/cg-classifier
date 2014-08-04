library("plyr")
library('dplyr')
library("VariantAnnotation")
library("tidyr")
library("reshape2")

process.genotypes <- function(g){
  
  clean <- function(x){
     if(is.numeric(x)) x[is.na(x)] = 0
    x
  }
  
  ## This creates a data.frame of relevant features
  ## for each genotype in the VCF file
  
  # this melts features in 1-dimensional matrices
  melt.1d = function(x) g[[x]] %>% 
    melt(varnames=c("VarID", "IID"), value.name=x)
  
  dfs.1d = lapply(c("GT", "FT", "GQ", "DP"), 
               melt.1d)
  
  # this melts features in 2-dimensional matrices
  melt.2d = function(x) g[[x]] %>%
    melt(varnames=c("VarID", "IID", "Allele")) %>%
    mutate(Allele=ifelse(Allele==1, paste(x,1,sep=""), paste(x,2,sep=""))) %>%
    dcast(VarID + IID ~ Allele)
  
  dfs.2d = lapply(c("AD", "HQ", "EHQ"),
                  melt.2d)
     
  ## bind all of the features together
  f.1d = do.call(cbind, lapply(dfs.1d, function(x) x[, 3, drop=F]))
  f.2d = do.call(cbind, lapply(dfs.2d, function(x) x[, 3:4, drop=F]))
  
  f = cbind(f.1d, f.2d)
  tmp = cbind(dfs.1d[[1]][ , 1:2], f)
  
  ## cleaning up for missing values
  tmp %>%
    lapply(clean) %>% 
    as.data.frame %>%
    mutate(FT=ifelse(FT %in% c("PASS", "VQLOW"), FT, "Other"))
}

process.sequences <- function(rd){
  
  ## This returns a data.frame of REF and ALT sequences
  ## and some basic comparisons between them
  
  ## this is messy, but the important part is a 
  ## series of mutation statements that add features
  ## to the data frame
  
  bases = elementMetadata(rd) %>%
    as.data.frame %>%
    dplyr::select(REF, ALT) %>%
    lapply(as.character) %>%
    as.data.frame(stringsAsFactors=F) %>% 
    mutate(REF.length=nchar(REF), ALT.length=nchar(ALT)) %>%
    mutate(Length.Diff=abs(REF.length-ALT.length)) %>%
    mutate(REF=ifelse(REF.length == 1, REF, "N")) %>% 
    mutate(ALT=ifelse(ALT.length == 1, ALT, "N")) %>%
    mutate(VarID=names(rd))
}

process.loci <- function(loci){
  
  collapse <- function(x){
    ## converts the compressed list to a normal list
    x = as.list(x)
    
    ## collapse loci with more than one entry
    ind = sapply(x, length) > 1
    x[ind] = lapply(x[ind], function(x) do.call(paste, as.list(x)))
    
    ## add entry if none
    ind = sapply(x, length) == 0
    x[ind] = c("-")
    
    unlist(x)
  }
  
  tmp = loci[c("CGA_FI", "CGA_RPT")] %>%
    lapply(collapse) %>%
    as.data.frame(stringsAsFactors=F) 
  
  tmp2 = tmp %>%
    mutate(IS.INTRON=grepl("INTRON", CGA_FI)) %>%
    mutate(IS.CDS=grepl("CDS", CGA_FI)) %>%
    mutate(IS.UTR=grepl("UTR", CGA_FI)) %>%
    mutate(IS.AT.RICH=grepl("AT_rich", CGA_RPT)) %>%
    mutate(IS.GC.RICH=grepl("GC_rich", CGA_RPT)) %>%
    mutate(IS.LOW.COMPLEXITY=grepl("Low_complexity", CGA_RPT)) %>%
    mutate(VarID=row.names(loci)) %>%
    dplyr::select(VarID, IS.INTRON, IS.CDS, IS.UTR, IS.AT.RICH, IS.GC.RICH, IS.LOW.COMPLEXITY)
}

vcfToDf <- function(v){
  
  ## Creates a data.frame with features for variant filter
  ## which can be use for either training data or testing data
  
  ## create data.frames for each type of data
  g = process.genotypes(geno(v))
  s = process.sequences(rowData(v))
  l = process.loci(info(v))
  
  ## merge data frames together
  merge(g, s) %>% merge(l)  
}