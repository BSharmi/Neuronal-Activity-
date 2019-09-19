rm(list=ls())

library(reshape)
library(ggplot2)
library(stringi)
library(WGCNA)
library(Rtsne)
library(RColorBrewer)
library(plyr)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

## enable multi thread for WGCBA
allowWGCNAThreads()

## for heatmap
#cols <- colorRampPalette(c('blue', 'yellow', 'red'))(50000)
#cols <- colorRampPalette(c('forestgreen', 'yellow', 'red'))(10)
#cols <- colorRampPalette(brewer.pal("YiOrRd"))(length(table.data))

## for complexheatmap. need to look at the figure and set these ranges
#col_fun = colorRamp2(c(0, 1, 4), c("forestgreen", "yellow", "red"))
col_fun = colorRamp2(c(0, 2), c("green", "red"))
col_fun(seq(-3, 3))

## setwd
setwd('/home/bsharmi6/NA_TF_project/RNAseq_dev/')

## read gene exp file
tpm.dat = read.table('TPM.all.txt', stringsAsFactors = F, sep = '\t', h=T)

## convert gsymb to upper
#gsymb = toupper(tpm.dat$Symbol)

## append/modify
tpm.dat$Symbol = toupper(tpm.dat$Symbol)

## merge for same genes
tpm.dat = ddply(tpm.dat, 'Symbol', numcolwise(mean, na.rm = T))

## select target genes of interest
#genelist = c('Egr1','Egr2', 'Klf14')
## retrieve rows from gene exp
#genes.exp = tpm.dat[grepl(paste(genelist, collapse="|"), x = tpm.dat$Symbol, ignore.case = F),]
## another way to retrieve genes exp
#genes.exp= tpm.dat[ tpm.dat$Symbol %in% genelist,]


## remove some columns
genes.exp = tpm.dat[, !grepl(paste(c('neurons', 'hind', 'mid', 'Fetal'), collapse="|"), colnames(tpm.dat),ignore.case = T)]
#rownames(genes.exp)=genes.exp$Symbol

#################################################################################### get TFs using CIS-BP ###########################################################
## read DBD file
mouse_DBD = read.delim(paste0('/home/bsharmi6/NA_TF_project/RNA-seq_all_TPM/CIS-BP-2.0.txt'), h = T, stringsAsFactors = F)
## keep relevant columns
mouse_DBD = mouse_DBD[, c('Name', 'Family')]
## convert to upper
mouse_DBD$Name = toupper(mouse_DBD$Name)
## split by comma
mouse_DBD = separate_rows(mouse_DBD, Family, sep = '\\,|\\/')

## get common in tpm semi join since we dont want duplicate genes and only TF is needed, not their domains
#tpm.annotated = semi_join(tpm.dat, mouse_DBD, by = c('Symbol' = 'mouseGene'))
genes.exp = semi_join(genes.exp, mouse_DBD, by = c('Symbol' = 'Name'))


##### mere replicates
.replicate.merged <- function(s) {
  gsymb = s[, 1]
  x = s[, 2:ncol(s)]
  colnames(x) <- gsub('^X', '', gsub('[.|_]CpG_ML.*', '', colnames(x)))
  colnames(x)[colnames(x) %in% "6wk_F"] <- "6wk"
  colnames(x)[colnames(x) %in% "6wk_neu_F"] <- "6wk NeuN+(F)"
  colnames(x)[colnames(x) %in% "7wk_neu"] <- "7wk NeuN+"
  colnames(x)[colnames(x) %in% "12mo_neu_F"] <- "12mo NeuN+(F)"
  colnames(x)[colnames(x) %in% "6wk_glia_F"] <- "6wk NeuN-(F)"
  colnames(x)[colnames(x) %in% "7wk_glia"] <- "7wk NeuN-"
  colnames(x)[colnames(x) %in% "12mo_glia_F"] <- "12mo NeuN-(F)"
  colnames(x)[colnames(x) %in% "7wk_glia"] <- "7wk NeuN-"
  colnames(x)[colnames(x) %in% "Excitatory.neuron"] <- "Excitatory Neu"
  colnames(x)[colnames(x) %in% "PV.neuron"] <- "PV Neu"
  colnames(x)[colnames(x) %in% "VIP.neuron"] <- "VIP Neu"
  colnames(x)[colnames(x) %in% "Neuron"] <- "E17.5 Neu"
  colnames(x) <-gsub("(m)(hind)|(m)(fore)|(m)(mid)", "\\2\\4\\6", colnames(x))
  
  if(!any(grepl(".rep[1-9]", colnames(x)))) {
    return(x)
  }
  
  y <- x[,!grepl(".rep[1-9]", colnames(x))]
  rep.names <- unique(gsub(".rep[1-9].*", "", colnames(x)[grepl(".rep[1-9]", colnames(x))]))
  
  z <- c()
  for(n in rep.names) {
    tmp <- x[,grepl(n, colnames(x))]
    if(is.null(dim(tmp))){
      z=cbind(z,tmp)
    }else{
      z <- cbind(z, rowSums(tmp, na.rm=T) / rowSums(!is.na(tmp)))
    }
  }
  colnames(z) <- rep.names
  
  w<-cbind(y, z)
  xlab.sorted = sort(colnames(w),index.return=T)
  xlab <- xlab.sorted$x; w <- w[,xlab.sorted$ix]
  lab.order = c(8:14, 7, 2, 4:6, 1, 3)  
  xlab <- xlab[lab.order]; w <- w[,lab.order]
  xlab = gsub('forebrain_', "", xlab)
  xlab = gsub('d0', 'P0', xlab)
  xlab = gsub('22mon', '22mo', xlab)
  ## add gsymb
  w<-cbind(gsymb, w)
  return(list(w,xlab))  
}


## call merge function
res<-.replicate.merged(genes.exp)
## get the gexp matrix
y_mat = res[[1]]
#y_mat = cbind.data.frame(Symbol = y_mat$gsymb, log2(1+y_mat[, 2:ncol(y_mat)]))
xlab = res[[2]]

## remove genes with 0's in atleast 1/3
y_mat = y_mat[apply(y_mat[, 2:ncol(y_mat)], 1, function(x) length(which(x==0))< 6),]
rownames(y_mat) = NULL

## divide by mean of each gene
#y_mat = cbind.data.frame(gsymb = y_mat[, 1], y_mat[,2:ncol(y_mat)]/rowMeans(y_mat[, 2:ncol(y_mat)]))

## normalize to log2
y_mat = cbind.data.frame(gsymb = y_mat[, 1], log2(y_mat[,2:ncol(y_mat)]+1))

########################################################################### kmeans ##################################
#cl <- kmeans(x = y_mat[, 2:ncol(y_mat)], centers = 12, iter.max = 50, nstart = 25)
#if (cl$ifault==4) {
# cl <- kmeans(x = y_mat[, 2:ncol(y_mat)], centers = 12, iter.max = 50, nstart = 25, algorithm="MacQueen")
#}
cl <- kmeans(x = y_mat[, 2:ncol(y_mat)], centers = 15, iter.max = 10000, nstart = 50, algorithm="MacQueen")
save(cl,file = paste0(getwd(), '/Kmeans/RNAseq_dev_Kmeans_result.RData'))

## check number of clusters
num_clust = unique(cl$cluster)

#### define output list
moduleLabels_list=vector('list', length(num_clust))
names(moduleLabels_list) = num_clust 

## assign dms to clusters
for(iclust in 1:length(num_clust)){
  clust_indx = which(cl$cluster %in% names(moduleLabels_list)[iclust])
  moduleLabels_list[[iclust]] = y_mat[clust_indx,]
}

## get gene names in each cluster
gene_list=vector('list', length(unique(cl$cluster)))
names(gene_list)=unique(cl$cluster)
for(iclust in 1:length(gene_list)){
  gene_list[[iclust]]=y_mat$gsymb[cl$cluster %in% names(gene_list)[iclust]]
}

## write target dms files 
sapply(names(gene_list), function (x) write.table(gene_list[[x]], file=paste0(getwd(), '/Kmeans/Kmeans_' ,x, "_genelist.txt"),append = F,quote = F,sep = '\t', row.names = F,col.names = F))

## box plots for  different time for each cluster of genes
for (iclust in 1:length(gene_list)){
  ##  retrieve genes exp
  genes.exp.iclust = y_mat[y_mat$gsymb %in% gene_list[[iclust]],]
  ## add expression for duplicate gene symbol
  genes.exp.iclust = ddply(genes.exp.iclust,"gsymb",numcolwise(mean))
  ## plot
  png(paste0(getwd(), '/Kmeans/Kmeans_gene_exp_dev_',names(gene_list)[[iclust]],'.png'), type = 'cairo', width = 650, height = 450); par(mai=c(3,1,1,1))
  pheatmap(genes.exp.iclust[, 2:ncol(genes.exp.iclust)], color = cols, cluster_col=F, show_rownames=F, clustering_method="ward.D2", 
           fontsize=18, cellwidth=35)
  dev.off()
}


####################################################################### WGCNA clustering ###########################################################################
net = blockwiseModules(t(y_mat[, 2:ncol(y_mat)]), power = 3, maxBlockSize = 50000, TOMType = "signed", minModuleSize = 200, reassignThreshold = 1e-6, mergeCutHeight = 0.15, numericLabels = TRUE, 
                       networkType = "signed", replaceMissingAdjacencies = TRUE, pamRespectsDendro = FALSE, loadTOM = FALSE, verbose = 3)
geneTree = net$dendrograms[[1]]
moduleLabels = net$colors
mergedColors = labels2colors(net$colors)
MEs = net$MEs
save(net, MEs, moduleLabels, geneTree,file = paste0(getwd(), '/WGCNA/RNAseq_dev_WGCNA_result.RData'))


moduleLabels = net$colors ## membership information of each DMS
mergedColors = labels2colors(net$colors)
unique(mergedColors)

## get modules
moduleLabels_list=vector('list', length(unique(moduleLabels)))
names(moduleLabels_list)=unique(moduleLabels)
for(i in 1:length(moduleLabels_list)){
  moduleLabels_list[[i]]=y_mat[moduleLabels %in% names(moduleLabels_list)[i],]
}

## get gene names in each cluster
gene_list=vector('list', length(unique(moduleLabels)))
names(gene_list)=unique(moduleLabels)
for(i in 1:length(gene_list)){
  gene_list[[i]]=y_mat$gsymb[moduleLabels %in% names(gene_list)[i]]
}
## write gene list
sapply(names(gene_list), function(x) write.table(gene_list[[x]], file = paste0(getwd(), '/WGCNA/', 'WGCNA_cluster_', x, '_genelist.txt'), append = F, quote = F, row.names = F, col.names = F))

## box plots for  different time for each cluster of genes
for (iclust in 1:length(gene_list)){
  ##  retrieve genes exp
  genes.exp.iclust = y_mat[y_mat$gsymb %in% gene_list[[iclust]],]
  ## add expression for duplicate gene symbol
  #genes.exp.iclust = ddply(genes.exp.iclust,"gsymb",numcolwise(sum))
  ## exclude gene esymb
  x <- genes.exp.iclust[, 2:ncol(genes.exp.iclust)]
  ### remove NA for plotting
  x <- x[rowSums(x<0, na.rm=T)<=0 & rowSums(!is.na(x))>=ncol(x),]
  ## rename
  colnames(x) = gsub('forebrain_','',colnames(x))
  ## plot
  png(paste0(getwd(), '/WGCNA/', 'WGCNA_gene_exp_dev_',names(gene_list)[[iclust]],'.png'), type = 'cairo', width = 450, height = 225); par(mai=c(3,1,1,1))
  #pheatmap(x, color = cols, cluster_col=F, show_rownames=F, clustering_method="ward.D2", symm=F,symkey=F, symbreaks=T, scale="none", fontsize=18, cellwidth=22)
  ht = Heatmap(x/rowMeans(x), gap = unit(2, "mm"),name = "log2(1+TPM)", col = col_fun,clustering_method_rows = "ward.D2",show_row_dend = F,cluster_columns = F,show_row_names = F,column_names_gp = gpar(fontsize = 20),
               heatmap_legend_param = list(legend_height = unit(6, "cm"),labels_gp = gpar(fontsize = 20),title_gp = gpar(fontsize = 20, fontface = "bold")))
  draw(ht)
  dev.off()
}

## assign cluster genes to three types based on expression
## set type
type = c('Postnatal', 'Mixed', 'Embryonic')

## combine all clustering together with a big heatmap
combined.cluster.dat = data.frame(matrix(NA,ncol = ncol(y_mat) +1))

## set col names
colnames(combined.cluster.dat) = c(colnames(y_mat), 'Type')

## fill up the matrix for different clusters
for(iclust in 1:length(gene_list)){
  #for(iclust in 1:length(moduleLabels_list)){
  ## get genes
  row_indx = which(y_mat$gsymb %in% gene_list[[iclust]])
  ## fill up matrix
  combined.cluster.dat = rbind.data.frame(combined.cluster.dat, cbind.data.frame(y_mat[row_indx,], Type = type[iclust]))
}

## remove first
combined.cluster.dat = combined.cluster.dat[-1,]

## set rownames to null
rownames(combined.cluster.dat) = NULL

## change to log
combined.cluster.dat = cbind.data.frame(Symbol = combined.cluster.dat[, 1], Type = combined.cluster.dat[,ncol(combined.cluster.dat)], combined.cluster.dat[, 2:(ncol(combined.cluster.dat)-1)])

## write
write.table(combined.cluster.dat, file = paste0(getwd(), '/WGCNA/', 'combined.cluster.dat.full.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = T)

## rename
colnames(combined.cluster.dat) <- gsub('forebrain_', '', colnames(combined.cluster.dat))

## melt data for ggplot
ggdat <- melt(combined.cluster.dat,id.vars = c("Symbol", "Type"))

## for finding top TFs in each time
max.tpm.symb = vector('list', ncol(combined.cluster.dat)-2)
names(max.tpm.symb) = colnames(combined.cluster.dat)[3:ncol(combined.cluster.dat)]
for(itime in colnames(combined.cluster.dat)){
  ## break if not relevant column
  if(itime %in% c('Type', 'Symbol'))
    next()
  else{
    max.tpm = sort(combined.cluster.dat[,itime], decreasing = T)[1:3]
    max.tpm.symb[[itime]] = combined.cluster.dat$Symbol[combined.cluster.dat[,itime] %in% max.tpm]
    #cat(max.tpm.symb)
    #cat('\n')
  }
}

## define colors
cbPalette=c("red","khaki", 'blue')
#cbPalette=c("forestgreen","antiquewhite", "antiquewhite1", "gold1", 'blue')

## plot
theme_set(theme_classic(base_size = 28)) 
png(paste0(getwd(), '/WGCNA/', 'gene_exp_dev_barplot.png'), type = 'cairo', width = 1000, height = 800); par(mai=c(3,1,1,1))
ggplot(ggdat, aes(fill=Type, y=value, x=variable)) + geom_jitter(aes(colour = Type)) +  xlab("") + ylab("TPM(Log scale)") +
  scale_color_manual(values=cbPalette) +theme_set(theme_classic(base_size = 18)) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()


