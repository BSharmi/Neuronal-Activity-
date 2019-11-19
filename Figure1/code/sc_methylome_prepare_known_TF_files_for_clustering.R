library(plyr)
library(dplyr)
library(data.table)
library(tidyverse)


############################################## selecting regions from TF motif from subset dms #######################################
rm(list=ls())

## define path
fpath = '/home/bsharmi6/NA_TF_project/scMethylome/Clustering_in_ETRM/subset_DMS/'

## list TF folders
TFlist <- list.files(fpath)

## read all dms all samples file
all.dms.samples.dat <- read.table(paste0('/home/bsharmi6/NA_TF_project/scMethylome/all.dms.all.samples.txt'), h=T, stringsAsFactors=F)

## select only 16 neurons
all.dms.samples.dat <- all.dms.samples.dat[, c(1:3,78:ncol(all.dms.samples.dat))]

## read all dms background
all.dms.back.dat <- read.table(paste0('/home/bsharmi6/NA_TF_project/scMethylome/all_DMS_known_TFs/all.dms_background.txt'), h=T, sep = '\t', stringsAsFactors=F)

## get regions for each TF and prepare files for clustering
for(iTF in 1:length(TFlist)){
	## read motif file
	TF.dat <- read.table(paste0(fpath, TFlist[iTF],'/', TFlist[iTF],'.motif.txt'), h=T, sep = '\t')
	## split first column
	TF.dat <- separate(TF.dat, FASTA.ID, c('chrom', 'start', 'end'), sep = ":|-") %>% select(chrom, start, end) %>% mutate_at(c("start", "end"), as.numeric) %>% unique()
	## get intersection of regions in all dms
	TF.meth.dat <- inner_join(TF.dat, all.dms.samples.dat, by=c('chrom', 'start', 'end'))
	## get background
	TF.back.dat <- all.dms.back.dat[sample(nrow(all.dms.back.dat), 2*nrow(TF.meth.dat)),]
	## write dms, meth and background for clustering
	write.table(TF.meth.dat[, 1:3], file = paste0(fpath, TFlist[iTF], '/', TFlist[iTF],'_dms.txt'),sep = '\t', col.names = T, row.names = F, quote = F, append = F)
	write.table(TF.back.dat, file = paste0(fpath, TFlist[iTF], '/', TFlist[iTF],'_background.txt'),sep = '\t', col.names = T, row.names = F, quote = F, append = F)
	write.table(TF.meth.dat[, 4:ncol(TF.meth.dat)], file = paste0(fpath, TFlist[iTF], '/', TFlist[iTF],'_methylation.matrix.txt'),sep = '\t', col.names = T, row.names = F, quote = F, append = F)
}

################################################# selecting regions from TF motif in all of dms ####################################
rm(list=ls())

## define path
fpath = '/home/bsharmi6/NA_TF_project/scMethylome/all_DMS_known_TFs/'

## list TF folders
TFlist <- list.files(fpath, pattern = 'minus', ignore.case=T)

## remove minus
TFlist <- gsub('_minus', '', TFlist)

## read all dms all samples file
all.dms.samples.dat <- read.table(paste0('/home/bsharmi6/NA_TF_project/scMethylome/all.dms.all.samples.txt'), h=T, stringsAsFactors=F)

## read all dms background
all.dms.back.dat <- read.table(paste0('/home/bsharmi6/NA_TF_project/scMethylome/all_DMS_known_TFs/all.dms_background.txt'), h=T, sep = '\t', stringsAsFactors=F)

## get regions for each TF and prepare files for clustering
for(iTF in 1:length(TFlist)){
	## read motif file
	TF.dat <- read.table(paste0(fpath, TFlist[iTF],'.motif.txt'), h=T, sep = '\t')
	## split first column
	TF.dat <- separate(TF.dat, FASTA.ID, c('chrom', 'start', 'end'), sep = ":|-") %>% select(chrom, start, end) %>% mutate_at(c("start", "end"), as.numeric) %>% unique()
	## get intersection of regions in all dms
	TF.meth.dat <- inner_join(TF.dat, all.dms.samples.dat, by=c('chrom', 'start', 'end'))
	## get background
	TF.back.dat <- all.dms.back.dat[sample(nrow(all.dms.back.dat), 2*nrow(TF.meth.dat)),]
	## write dms, meth and background for clustering
	write.table(TF.meth.dat[, 1:3], file = paste0(fpath, TFlist[iTF], '/', TFlist[iTF],'_dms.txt'),sep = '\t', col.names = T, row.names = F, quote = F, append = F)
	write.table(TF.back.dat, file = paste0(fpath, TFlist[iTF], '/', TFlist[iTF],'_background.txt'),sep = '\t', col.names = T, row.names = F, quote = F, append = F)
	write.table(TF.meth.dat[, 4:ncol(TF.meth.dat)], file = paste0(fpath, TFlist[iTF], '/', TFlist[iTF],'_methylation.matrix.txt'),sep = '\t', col.names = T, row.names = F, quote = F, append = F)
}