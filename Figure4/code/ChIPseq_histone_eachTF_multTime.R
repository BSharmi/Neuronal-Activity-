rm(list=ls())
library(ggplot2)
library(foreach)
library(doParallel)
library(dplyr)
library(tidyr)

## parallel processing
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#cl <- makeCluster(20)
registerDoParallel(cl)


## define path
fpath = '/home/bsharmi6/NA_TF_project/HistoneMark/Forebrain/'
## define outputpath if not already present
if(!dir.exists(paste0(fpath, 'results'))){
  dir.create(file.path(paste0(fpath, 'results'))) 
}
## define time points
timepoints = dir(paste0(fpath, 'data'))
## delete 10.5 because not all histone are present
timepoints = timepoints[-2]

## define TFs
TFname = 'MEF2A'

## funciton to calculate overlap
get_overlap = function(dat_f, hist_f){
  res<- NA
  dat.start = as.numeric(dat_f["start"]); dat.end = as.numeric(dat_f["end"]); dat.chr = as.character(dat_f["chrom"]); dat.PeakID = as.character(unname(as.matrix(dat_f["PeakID"])))
  hist_peak_indx = as.character(unname(as.matrix(dat_f[4])))
  hist.start = as.numeric(hist_f[hist_f[,4] %in% hist_peak_indx, 2])
  hist.end = as.numeric(hist_f[hist_f[,4] %in% hist_peak_indx, 3])
  hist.chr = hist_f[hist_f[,4] %in% hist_peak_indx, 1]
  ## check if null overlap then set to 0
  if(hist_peak_indx == "Peak_0"){
    res = list(overlap=0, peak = data.frame(chrom = dat.chr, start = dat.start, end = dat.end, PeakID = dat.PeakID))
  }else if(hist.chr == dat.chr){
    ## new logic
    if(dat.start <  hist.end & dat.end > hist.start){
      res = list(overlap = as.numeric(dat.end - hist.start), peak = data.frame(chrom = dat.chr, start = dat.start, end = dat.end, PeakID = dat.PeakID))
    }
  }else {
    res <- NA
  }  
  return(res)
}

## for the TF
TF.list.timepoint = vector('list', length(timepoints))
#names(TF.list) = sapply(timepoints, function(x) strsplit(x, '\\.')[[1]][1])
names(TF.list.timepoint) = timepoints
for (itime in 1:length(timepoints)){
  ## for each timepoint
  dat = read.delim(paste0(fpath, '/data/', timepoints[itime], '/', TFname, '/', TFname, '.bed.annotated.txt'), h=T,stringsAsFactors = F)
  ## remove redundant rows and columms
  dat <- dat[,!colnames(dat) %in% 'X']
  dat <- dat[!dat$chrom %in% 'chrom',]
  ## get TF total length
  dat.length = sum(dat$end-dat$start)
  ## get number of histone
  hist_names = colnames(dat)[-c(1:3)]
  ## add peak ID column
  dat$PeakID = paste0('Peak', seq(1, nrow(dat)))
  
  
  #### parallel 
  hist_overlap_each_time <-foreach (ihist = 1:length(hist_names),.combine = 'list',.multicombine = TRUE, .maxcombine = length(hist_names),.errorhandling = "pass", .verbose=TRUE, .final = function(x) setNames(x, hist_names), .packages = c('tidyr')) %dopar% {
    ## old code discarding null overlap peaks
    #hist.indx = which(dat[,colnames(dat) %in% hist_names[ihist]] != 'null')
    #tmp = cbind.data.frame(dat[hist.indx, 1:3], dat[hist.indx,colnames(dat) %in% hist_names[ihist]])
    ## get null rows
    hist.nullindx = which(dat[,colnames(dat) %in% hist_names[ihist]] == 'null')
    ## set null overlap to an arbitrary peak
    dat[hist.nullindx, colnames(dat) %in% hist_names[ihist]] = 'Peak_0'
    tmp = cbind.data.frame(dat[, 1:3], HistPeak = dat[,colnames(dat) %in% hist_names[ihist]], PeakID = dat$PeakID)
    rownames(tmp) = NULL
    ## split multiple using tidyr shortest code option 3
    merged_dat = tmp %>% separate_rows(HistPeak, sep = ';')
    ## split multiple using tidyr
    ## split on ';' for gene names and duplicate the rows option 1
    # x = tmp %>%
    #   mutate(HistPeak = strsplit(as.character(HistPeak), ";")) %>%
    #   unnest() %>%
    #   filter(HistPeak != "") 
    ## split multiple using my own code option 2
    # dup.indx = grep('\\;',tmp$HistPeak)
    # ## check if any duplicates
    # if(length(dup.indx)>0){
    #   tmp.single = tmp[-dup.indx,]
    #   tmp.mult = tmp[dup.indx,]
    #   #tmp.mult.list <- split(tmp.mult, seq(nrow(tmp.mult)))
    #   ## split multiple peaks
    #   dup.peaks = sapply(tmp.mult[,4], function(x) strsplit(as.character(x), '\\;'))
    #   tmp.mult = tmp.mult[rep(row.names(tmp.mult), sapply(dup.peaks, length)),]
    #   tmp.mult[,4] = unlist(dup.peaks)
    #   merged_dat = rbind.data.frame(tmp.single, tmp.mult)
    #   rownames(merged_dat) = NULL
    # }else{
    #   merged_dat = tmp
    # }
    hist.dat = read.delim(paste0(fpath, '/data/' ,timepoints[itime], '/All_histones/', hist_names[ihist], '.bed'), h=F, stringsAsFactors = F)
    if(nrow(merged_dat)>0){
      ## call function
      apply(merged_dat, 1, function(x,y) get_overlap(x,hist.dat))
    }else 0   
  }
  
  ## remove NA or 0. these are with length 1
  hist_overlap_each_time <- lapply(hist_overlap_each_time, function(y) y[sapply(y, function(x) length(x) >1)])
  ## split overlap and peak
  hist.overlap = sapply(hist_overlap_each_time, function(x) sapply(x, '[[', 1))
  TF.peak = sapply(hist_overlap_each_time, function(x) sapply(x, '[[', 2, simplify = F))
  ## delete zero overlaps. NO DO NOT DELETE AS WE WANT TO KEEP THE 0 OVERLAP PEAKS. HOWEVER THE CODE IS GOOD FOR LATER REFERENCE
  #TF.peak = mapply(function(x,y) x[y], TF.peak, lapply(hist.overlap, function(x) x!=0))
  #hist.overlap = lapply(hist.overlap, function(x) {x[x!=0]})
  ## add to timepoint list
  TF.list.timepoint[[itime]] = list(TF.peak, hist.overlap)
}

## save
save(TF.list.timepoint, file = paste0(fpath, '/results/' ,TFname, '_histone_enrichment.RData'))

