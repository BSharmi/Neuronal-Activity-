library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)
library(stringr)

rm(list=ls())

## define path 
fpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/'

## define output path
outpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/Enrichment_motifs_in_neurons/TRM/'

## list TF folders
dms.list <- list.files(fpath, pattern = '^[m*]')

## for each dms
for(idms in 1:length(dms.list)){
  ## get ETRM folders
  TRMs <- gsub('_minus', '', list.files(paste0(fpath, dms.list[idms]), pattern = 'minus'))
  ## read motif files
  for(iTRM in 1:length(TRMs)){
    ## read
    tfmotif.dat <- read.table(paste0(fpath, dms.list[idms], '/',TRMs[iTRM], '.motif.txt'), h=T, stringsAsFactors=F, sep = '\t')
    ## split
    tfmotif.dat <- separate(tfmotif.dat, FASTA.ID, c('chrom', 'start', 'end'), sep = ":|-") %>% select(chrom, start, end) %>% mutate_at(c("start", "end"), as.numeric) %>% unique()
    ## get output TRM file name
    out_trm_name <- read.delim(paste0(fpath, dms.list[idms], '/',TRMs[iTRM], '/knownResults.txt'), h=T, stringsAsFactors=F, sep = '\t' )[1,1]
    ## delete extra characters
    out_trm_name <- strsplit(out_trm_name,'\\(|\\:|\\.|\\-')[[1]][1]
    #cat(out_trm_name)
    #cat('\n')
    ## write and append
    write.table(tfmotif.dat, file = paste0(outpath,out_trm_name,'.trm.txt'), sep = '\t', col.names = !file.exists(paste0(outpath,out_trm_name,'.trm.txt')), row.names = F, quote = F, append = T)
    #if(file.exists(paste0(outpath,out_trm_name,'.trm.txt'))){
    #	write.table(tfmotif.dat, file = paste0(outpath,out_trm_name,'.trm.txt'), sep = '\t', col.names = F, row.names = F, quote = F, append = T)
    #}else{
    #	write.table(tfmotif.dat, file = paste0(outpath,out_trm_name,'.trm.txt'), sep = '\t', col.names = T, row.names = F, quote = F, append = F)
    #}
    
  }
}