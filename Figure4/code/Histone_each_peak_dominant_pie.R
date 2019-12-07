rm(list=ls())

library(plyr)
#library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)

## set wd
setwd('/home/bsharmi6/NA_TF_project/HistoneMark/Forebrain/results/')

## give TF name 
TF = 'Egr1'

## load data
load(paste0(TF, '_histone_enrichment.RData'))

## split peak and overlap into two lists
TF.peak = sapply(TF.list.timepoint, '[' , 1)
hist.overlap = sapply(TF.list.timepoint, '[' , 2)

## rbind TF.peak
TF.peak = lapply(TF.peak, function(x) sapply(x, function(y) do.call('rbind', y), simplify = F))

## annoate TF peaks with the maximum overlapped histone mark for each time point
timepoints = names(TF.peak)
## output list
peak.hist.overlap_res = vector('list', length(timepoints))
names(peak.hist.overlap_res) = timepoints
## run loop for each time
for(itime in timepoints){
  ## add histone overlap column to TF peak df
  tmp = mapply(function(x,y) cbind.data.frame(x,y), TF.peak[[itime]], hist.overlap[[itime]], SIMPLIFY = F)
  ## merge unique columns by adding overlap
  annotated.TF.dat = lapply(tmp, function(x) ddply(x, c("chrom", "start", "end", "PeakID"),numcolwise(sum)))
  ## merge annotated.TF.dat from a list to a single data frame by appending the overlap columns after the chrom, start, end, PeakID columns through full join
  merged.TF.dat = annotated.TF.dat %>% reduce(full_join, by = c("chrom", "start", "end", "PeakID"))
  ## add column nanses
  colnames(merged.TF.dat)[5:ncol(merged.TF.dat)] = names(tmp)
  ## assign None to no overlap peaks
  merged.TF.dat$hist_anno = apply(merged.TF.dat[, 5:ncol(merged.TF.dat)], 1, function(x) {ifelse(all(x ==0), 'None', list(names(x)[which(x==max(x))]))})
  peak.hist.overlap_res[[itime]] = table(unlist(merged.TF.dat$hist_anno))
}
save(peak.hist.overlap_res, file = paste0(getwd(),'/',TF, '_histone_overlap_fraction.RData'))
## rename d0 to P0
#names(peak.hist.overlap) = gsub('d0', 'P0', names(peak.hist.overlap), ignore.case = T)
## arrange
peak.hist.overlap_res= peak.hist.overlap_res[c(2:length(peak.hist.overlap_res), 1)]
## convert to fraction
peak.hist.overlap = lapply(peak.hist.overlap_res, function(x) as.double(format(x/nrow(merged.TF.dat) *100, nsmall = 1,digits = 1, trim = T)))
## ggplot
# ggdat = data.frame(Value = unlist(peak.hist.overlap), 
#                    Histone = unlist(sapply(peak.hist.overlap_res, function(x) strsplit(names(unlist(x)), '\\.H'))),
#                    Time = rep(names(peak.hist.overlap), length(peak.hist.overlap[[1]])))

# ggplot(ggdat, aes(x=Time, y=Value, group=Histone)) +
#   geom_line(aes(color=Histone))+
#   geom_point(aes(color=Histone)) + scale_color_brewer(palette="Dark2") +theme_classic()
# ## plot pie
# par(mfrow=c(2,4))
# for(i in 1:length(peak.hist.overlap_res)){
#   pielabels = paste0(names(peak.hist.overlap_res[[i]]), '=', peak.hist.overlap[[i]], '%' )
#   pie(peak.hist.overlap[[i]], labels=NA, clockwise=TRUE, col = brewer.pal(length(peak.hist.overlap[[i]]),"Set1"), border="white",
#       radius = 0.9, cex=0.8,main = names(peak.hist.overlap_res[i])) 
#   legend('bottomright',legend=pielabels, bty="n", fill=brewer.pal(length(peak.hist.overlap[[i]]),"Set1"))
# }
# legend(x = "right",inset = 1,
#        legend = pielabels, bty="n", fill=brewer.pal(length(peak.hist.overlap[[i]]),"Set1"),
#        lwd=5, cex=.5, horiz = TRUE)
#legend(1.6,1,legend=pielabels, bty="n", fill=brewer.pal(length(peak.hist.overlap[[i]]),"Set1"))

#par(mfrow=c(1,1))
# for(i in 1:length(peak.hist.overlap_res)){
#   png(paste0(getwd(),'/',TF, '_',names(peak.hist.overlap_res[i]),'.png'), type="cairo", width = 741, height = 545)
#   pielabels = paste0(names(peak.hist.overlap_res[[i]]), '=', peak.hist.overlap[[i]], '%' )
#   pie(peak.hist.overlap[[i]], labels=NA, clockwise=TRUE, col = brewer.pal(length(peak.hist.overlap[[i]]),"Set1"), border="white", radius = 0.5) 
#   legend(0.5,0.5,legend=pielabels, bty="n", fill=brewer.pal(length(peak.hist.overlap[[i]]),"Set1"),cex = 1.4) 
#   mtext(names(peak.hist.overlap_res[i]), side = 3,line=-3.5, cex=1.5)
#   dev.off() 
# }
