## TF motif enrichment
rm(list=ls())

library(pheatmap)
library(dplyr)
library(stringr)
library(Rtsne)
library(gplots)

## set annotation type
#ann.type = 'stable'

## define fpath 
#fpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMs/Enrichment_motifs_in_neurons/all_motif_locations/DMS/'
fpath = paste0('/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/Enrichment_motifs_in_neurons/DMS/')

## read annotated table. it is annotated with clusters and TF peaks
tmp = read.delim(paste0(fpath, list.files(fpath, pattern = '*annotated.txt')), header = T)

## create a column with cluster annotation
tmp <- rename(tmp, cluster_annotated = DMS)

# ## get trm column names
colnames(tmp) = gsub('.trm', '', colnames(tmp), ignore.case=T)

## set excitatory
tmp$cluster_annotated = gsub("mDL-2|mDL-1|mDL-3|mL23|mL4|mL6-2|mL6-1|mL5-1|mL5-2|mIn-1", "excitatory", tmp$cluster_annotated)
tmp$cluster_annotated = gsub("mPv|mSst-1|mNdnf-2|mNdnf-1|mSst-2|mVip", "inhibitory", tmp$cluster_annotated)
# ## set names for each
# for(iTRM in 1:length(trm_cols)){
# 	tmp$cluster_annotated[! tmp[trm_cols[iTRM]] == 0] = paste0(str_to_title(strsplit(trm_cols[iTRM], '\\.')[[1]][1]), '_', 'ETRM')
# }

# ## remove no TMRs
# tmp = tmp[! tmp$cluster_annotated %in% 'Non-TRMs',]

## delete trm columns
#tmp = tmp[,c(1:3,64,4:55)]
#tmp <- tmp[, -grep('.trm', colnames(tmp), ignore.case = FALSE)]
#tmp = tmp %>% select(1:3, ncol(tmp), everything())


## rename 1:3
colnames(tmp)[1:3] = c('chrom', 'start', 'end')
## remove .motifs
#colnames(tmp) <- gsub('.motifs.HOMER', '', colnames(tmp))
## remove end space
colnames(tmp) <- gsub('\\.$', '', colnames(tmp))
## capitalize
colnames(tmp) <- str_to_title(colnames(tmp))
## add ETRM extenstion to name
#colnames(tmp)[5:ncol(tmp)] <- paste0(colnames(tmp)[5:ncol(tmp)], '_ETRM')

#tmp$cluster_annotated[! tmp$cluster_5_dms %in% 0] = 'cluster_5'
## remove cluster 0
#tmp = tmp[! tmp$cluster_annotated %in% '0',]

#tmp = tmp[,c(1:3,52,4:48)]
# if(TF =='Tet1'){
# 	## order for Tet1
# 	tmp = tmp[,c(1:3,60,4:55)]
# }else{
# 	## order for Tet1
# 	tmp = tmp[,c(1:3,59,4:55)]
# }

## Fishers test
x=tmp
## factor by cluster dms
dat <- c(); total <- c()
for(f in levels(factor(x[,4]))) {
  t <- x[grepl(f, x[,4]),]
  dat <- rbind(dat, colSums(x[grepl(f, x[,4]),5:ncol(x)]))
  total <- rbind(total, sum(grepl(f, x[,4])))
  rm(list=c("t"))
} 
rownames(dat) <- levels(factor(x[,4]))
rownames(total) <- levels(factor(x[,4]))


dat.pval <- c(); dat.est <- c()
for(i in 1:ncol(dat)) for(j in 1:nrow(dat)) {
  tab <- matrix(c(dat[j,i], total[j], sum(dat[,i]), sum(total)), ncol=2) 
  f.et <- fisher.test(tab)
  dat.est <- c(dat.est, f.et$estimate)
  dat.pval <- c(dat.pval, f.et$p.value)
}

dat.est <- matrix(dat.est, ncol=ncol(dat), dimnames=list(rownames(dat), colnames(dat)))
dat.pval <- matrix(dat.pval, ncol=ncol(dat), dimnames=list(rownames(dat), colnames(dat)))
#df.dat <- dat.est[,setdiff(colnames(dat.est), c("SRF", "CREB", "Miz1", "Con_enhancers", "Dec_enhancers", "Inc_enhancers", "H3K27ac.KCl","Enhancer"))]
df.dat <- dat.est

## dirty approach to shorten NF1 column name
#colnames(df.dat) <- gsub('.Halfsite', '', colnames(df.dat))
## change cluster names to bring egr1 to front
#if(TF == 'Tet1'){
#rownames(df.dat) = c("cluster_3", "cluster_1", "cluster_2", "cluster_4")
#}else{
#	rownames(df.dat) = c("cluster_2", "cluster_1", "cluster_4", "cluster_3")
#}

## order rows since not clustering on rows
#df.dat = df.dat[order(rownames(df.dat)),]
if(nrow(df.dat) > 2){
  cols <- colorRampPalette(c('green', 'white', 'red'))(100)
  tiff(paste0(fpath, "sixteen_neuron.TF.DMS.corr.tiff"), type="cairo", width = 2000, height = 1000)
  #pheatmap(cor(df.dat), color = cols, show_rownames=T, clustering_method="average", fontsize=18)
  pheatmap(cor(df.dat), color = cols, show_rownames=T, cluster_cols = T, cluster_rows = F, clustering_method="ward.D2", fontsize=22)
  dev.off()
  
  ## neuron not needed if we dont see clusters showing difference
  #excitatory.neurons <- c("cluster_1")
  #inhibitory.neurons <- c("cluster_2")
  #neuron <- rep(NA, nrow(dat.pval))
  #neuron[rownames(dat.est) %in% excitatory.neurons] <- "Excitatory"
  #neuron[rownames(dat.est) %in% inhibitory.neurons] <- "Inhibitory"
  
  #anno <- data.frame(Neuron = factor(neuron))
  #anno_color <- list(Neuron = c("orange1", "cyan"))
  #names(anno_color$Neuron) <- levels(anno$Neuron)
  #rownames(anno) <- rownames(dat.est)
  
  rate <- (max(df.dat)-1)/(max(df.dat)-min(df.dat))
  cols <- c(colorRampPalette(c('green', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'red'))(round(rate*100)))
  ## tiff
  tiff(paste0(fpath,"sixteen_neuron.TF.FC.tiff"), type="cairo", width = 2000, height = 800)
  par(mai=c(1,1,1,1))
  #pheatmap(df.dat, color = cols, annotation_row=anno, annotation_colors=anno_color, show_rownames=T, clustering_method="ward.D2", fontsize=18, cellwidth=22, , cellheight=20)
  pheatmap(df.dat, color = cols, show_rownames=T, clustering_method="ward.D2", cluster_rows = FALSE, cluster_cols = T, fontsize=18, cellwidth=22, cellheight=20)
  dev.off()
  
  ## png
  rate <- (max(log2(df.dat+0.2))-1)/(max(log2(df.dat+0.2))-min(log2(df.dat+0.2)))
  cols <- c(colorRampPalette(c('green', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'red'))(round(rate*100)))
  png(paste0(fpath, "sixteen_neuron.TF.FC",".png"), type="cairo", width = 800, height = 600);#par(mar=c(7,8,7,10))#
  #par(mai=c(3,1,1,1))
  #pheatmap(df.dat, color = cols, annotation_row=anno, annotation_colors=anno_color, show_rownames=T, clustering_method="ward.D2", fontsize=18, cellwidth=22, , cellheight=20)
  pheatmap(log2(df.dat+0.2), color = cols, show_rownames=T, clustering_method="ward.D2", cluster_rows = T, cluster_cols = T, fontsize=26, cellwidth=22, cellheight=24)
  dev.off()
  
  ## png with colored rownames
  rate <- (max(log2(df.dat+0.2))-1)/(max(log2(df.dat+0.2))-min(log2(df.dat+0.2)))
  cols <- c(colorRampPalette(c('green', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'firebrick1'))(round(rate*100)))
  ann_colors <- c(rep("Blue", 10), rep("orange3", 6))
  png(paste0(fpath, "sixteen_neuron.TF.FC",".png"), type="cairo", width = 1100, height = 900);par(oma=c(8,1,1,10))#par(mai=c(5,0.2,0.2,5))#
  heatmap.2(as.matrix(log2(df.dat+0.2)), scale = "none", col = cols, hclustfun = function(x) hclust(x, method="ward.D2"), colRow=ann_colors, cexCol=3, cexRow=3, trace = "none", density.info = "none", key = F,dendrogram='none', Rowv=TRUE, Colv=TRUE,labRow = rownames(df.dat),
            sepwidth=c(0.5,0.08), offsetRow = 0.01, offsetCol = 0.01, pty = 'm',adjCol = c(1.0,0.5), adjRow = c(0.05, 0.3)) 
  legend(y=0.29, x=0.05, xpd=TRUE, legend = seq(-2,2),col = colorRampPalette(c('green', 'white','red'))(5), lty= 1,lwd = 5, cex=1.7)
  dev.off()
}else{
  cols <- colorRampPalette(c('green', 'white', 'red'))(100)
  tiff(paste0(fpath, "sixteen_neuron.TF.DMS.corr.grouped.tiff"), type="cairo", width = 2000, height = 1000)
  #pheatmap(cor(df.dat), color = cols, show_rownames=T, clustering_method="average", fontsize=18)
  pheatmap(cor(df.dat), color = cols, show_rownames=T, cluster_cols = T, cluster_rows = F, clustering_method="ward.D2", fontsize=22)
  dev.off()
  
  ## neuron not needed if we dont see clusters showing difference
  #excitatory.neurons <- c("cluster_1")
  #inhibitory.neurons <- c("cluster_2")
  #neuron <- rep(NA, nrow(dat.pval))
  #neuron[rownames(dat.est) %in% excitatory.neurons] <- "Excitatory"
  #neuron[rownames(dat.est) %in% inhibitory.neurons] <- "Inhibitory"
  
  #anno <- data.frame(Neuron = factor(neuron))
  #anno_color <- list(Neuron = c("orange1", "cyan"))
  #names(anno_color$Neuron) <- levels(anno$Neuron)
  #rownames(anno) <- rownames(dat.est)
  
  rate <- (max(df.dat)-1)/(max(df.dat)-min(df.dat))
  cols <- c(colorRampPalette(c('steelblue', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'red'))(round(rate*100)))
  ## tiff
  tiff(paste0(fpath,"sixteen_neuron.TF.FC.grouped.tiff"), type="cairo", width = 2000, height = 800)
  par(mai=c(1,1,1,1))
  #pheatmap(df.dat, color = cols, annotation_row=anno, annotation_colors=anno_color, show_rownames=T, clustering_method="ward.D2", fontsize=18, cellwidth=22, , cellheight=20)
  pheatmap(df.dat, color = cols, show_rownames=T, clustering_method="ward.D2", cluster_rows = FALSE, fontsize=18, cellwidth=22, cellheight=20)
  dev.off()
  
  ## arrange manually if cluster_cols = F; BUT NOT A GOOD APPROACH
  df.dat <- df.dat[,order(df.dat[2,], decreasing = T)] ## ordering by inhibitory; can also be by excitatory
  #df.dat <- df.dat[, c("Tcf21", "Mafa", "Ap4", "NKX6.1", "Rorgt", "Lhx3", "Sox3", "Egr1", "Atoh1", "Junb", "Tbr1", "Neurog2", "Rfx6", "Atf3", "Nur77", "Fra1", "Tgif2", "Mef2b", "Nf1", "Oct4")]
  ## png
  #rate <- (max(log2(df.dat+0.1))-1)/(max(log2(df.dat+0.1))-min(log2(df.dat+0.1)))
  #cols <- c(colorRampPalette(c('steelblue', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'red'))(round(rate*100)))
  png(paste0(fpath, "sixteen_neuron.TF.FC.grouped.png"), type="cairo", width = 800, height = 600);#par(mar=c(7,8,7,10))#
  #par(mai=c(3,1,1,1))
  #pheatmap(df.dat, color = cols, annotation_row=anno, annotation_colors=anno_color, show_rownames=T, clustering_method="ward.D2", fontsize=18, cellwidth=22, , cellheight=20)
  pheatmap(log2(df.dat+0.1), color = cols, show_rownames=T, clustering_method="ward.D2", cluster_rows = F, cluster_cols = F, fontsize=26, cellwidth=22, cellheight=24)
  dev.off()
  
  ## png with colored rownames
  cols <- c(colorRampPalette(c('dodgerblue2', 'white'))(100-round(rate*100)), colorRampPalette(c('white', 'firebrick1'))(round(rate*100)))
  png(paste0(fpath, "sixteen_neuron.TF.FC.grouped.png"), type="cairo", width = 800, height = 600);par(oma=c(8,1,1,10))#par(mai=c(5,0.2,0.2,5))#
  heatmap.2(as.matrix(log2(df.dat+0.1)), scale = "none", col = cols, cexCol=3, cexRow=3, trace = "none", density.info = "none", key = F,dendrogram='none', Rowv=F, Colv=F,labRow = rownames(df.dat),
            offsetRow = 0.01, offsetCol = 0.01, pty = 'm',adjCol = c(1.0,0.5), adjRow = c(0.05, 0.3), margins = c(25,2)) 
  dev.off()
}

