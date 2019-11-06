rm(list=ls())
library(pheatmap)
library(RColorBrewer)
library(WGCNA)
library(plyr)
library(dplyr)
library(tidyr)
#library(data.table)

## for heatmap
cols <- colorRampPalette(c('blue', 'yellow'))(100)

## define path
fpath='/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/ChIP_seq_analysis/' 

## define TF
TF = 'Egr1'

## define sclist
sclist = list.files(paste0(fpath, TF), pattern = '^[m*]', ignore.case=F)


##### function for plotting
.plot.profile <- function(t,iclust_DMS,xlab,clusterrows_func=T,cex=1) {
  ### remove NA for plotting
  x <- t[rowSums(t<0, na.rm=T)<=0 & rowSums(!is.na(t))>=ncol(t),]
  ############### heatmap
  #pdf(paste0('DMS_clust',iclust_DMS,'_heatmap.pdf'), width = 15, height = 13); par(mai=c(3,1,1,1))
  png(paste0('DMS_clust_',iclust_DMS,'_heatmap_customized.png'),type="cairo", width = 800, height = 600); par(mai=c(3,1,1,1))
  pheatmap(x, color = cols, cluster_col=T, cluster_rows = clusterrows_func, show_rownames=F, clustering_method="ward.D2",fontsize=35, cellwidth=35)
  dev.off()
  ############## methylation profile
  apply(x, 2, mean, na.rm=T) -> avg
  apply(x, 2, sd, na.rm=T) -> sdev
  upper <- approx(1:ncol(x), avg+sdev, n=100, method = "linear")
  lower <- approx(1:ncol(x), avg-sdev, n=100, method = "linear")
  #pdf(paste0('DMS_clust',iclust,'.pdf'),width = 15, height = 13)
  png(paste0('DMS_clust_',iclust_DMS,'_profile_customized.png'),type="cairo", width = 1040, height = 580)
  par(mai=c(2.5,1.5,0.82,0.42))
  plot(1:ncol(x), avg, type='n', ylim=c(min(lower$y), max(upper$y)), xlab='', ylab='Methylation level', xaxt='n', yaxt='s', cex.axis=cex*1.8, cex.lab=cex*1.8, cex=cex*1.8)
  polygon(c(rev(lower$x), upper$x), c(rev(lower$y), upper$y), col = 'grey80', border = NA)
  ## color points by cluster color
  #points(1:ncol(x), avg, pch=16, col=coltype[iclust_DMS], cex=1.6*cex)
  ## color black for all clusters
  points(1:ncol(x), avg, pch=16, col='black', cex=1.6*cex)
  lines(1:ncol(x), avg, col='black', lty='dashed', lwd=1.8*cex)
  axis(1, at=seq_along(1:ncol(x)), label=xlab, las=2, cex.axis=1.8, cex.lab=cex*1.8, cex=cex*1.8)
  dev.off()
}

## plot for all TFS
for (isc in 1:length(sclist)){
  #setwd('/home/bsharmi6/Methyome_data_analysis/TET1KO_TET2KO/TET1KO/DMS_6week_analysis/DMR_up_down_background/Kmeans_WGCNA/Kmeans_cluster_3/')
  setwd(paste0(fpath,TF,'/',sclist[isc]))
  ## read meth matrix
  y <- read.table(paste0(fpath,TF,'/',sclist[isc], '/',sclist[isc], '_methylation.matrix.txt'), h=T, sep = '\t')
  ## change column names
  colnames(y) <- gsub('_CpG_ML', '', colnames(y))
  ## read.dms
  tmp = y
  xlab = colnames(y)
  .plot.profile(tmp,sclist[isc],xlab, clusterrows=T)
}
