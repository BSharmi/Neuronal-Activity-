## selecting cell type specific brain RNA-seq data

## add libraries
library(readxl)
library(plyr)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape)
library(homologene)
library(kmed)
library(biclust)

rm(list=ls())

## for heatmap
cols <- colorRampPalette(c('forestgreen', 'yellow', 'red'))(10)

## read annotation file
mm10.gtf = read.table('/home/bsharmi6/Annotation/mm10.genes.gtf.txt', h=T, stringsAsFactors = F)

## define path 
fpath = '/groups/ECBL/Xiaoran_codes_testing/RNA_seq_cell_type/'

## list directories
filelist <- list.files(fpath, pattern = 'SRR')

## read the first file and get genes
tmp = read.table(dir(paste0(fpath,filelist[1]), full.names=T, pattern = 'genes.results$'), h=T, stringsAsFactors = F)
## get sample name
samplename = tail(strsplit(dir(paste0(fpath,filelist[1]), full.names=T, pattern = 'genes.results$'), split = '\\/')[[1]], n=1)
## get gene names
tpm.dat = inner_join(tmp, mm10.gtf, by = c('gene_id' = 'gene_id'))[, c("gene_symbol", "TPM")]
## colname change
colnames(tpm.dat)[2] = samplename


## read gene file
for(ifile in 2:length(filelist)){
	## read a file
	tmp = read.table(dir(paste0(fpath,filelist[ifile]), full.names=T, pattern = 'genes.results$'), h=T, stringsAsFactors = F)
	## get sample name
	samplename = tail(strsplit(dir(paste0(fpath,filelist[ifile]), full.names=T, pattern = 'genes.results$'), split = '\\/')[[1]], n=1)
	## get common genes with gene table
	common_dat = inner_join(tmp, mm10.gtf, by = c('gene_id' = 'gene_id'))[, c("gene_symbol", "TPM")]
	## check if all same
	if(all(tpm.dat$gene_symbol ==  common_dat$gene_symbol)){
		## append
		tpm.dat = cbind.data.frame(tpm.dat, common_dat[common_dat$gene_symbol %in% tpm.dat$gene_symbol,2])
		## rename column
		colnames(tpm.dat)[ifile +1] = samplename
	}else{
		next()
	}
}

## shorten column names
colnames(tpm.dat) = gsub(".genes.results", "", colnames(tpm.dat))

## append/modify
tpm.dat$gene_symbol = toupper(tpm.dat$gene_symbol)

## remove WC
tpm.dat = tpm.dat[, !grepl('WC', colnames(tpm.dat))]

## merge for same genes
tpm.dat = ddply(tpm.dat, 'gene_symbol', numcolwise(mean, na.rm = T))

########################################################################################## using homologene package ##########################################
# ## read TF db file
# human.tf.dat = read_excel(paste0(fpath, 'The_Human_Transcription_Factors_full_database.xlsx'), sheet = 1)

# ## change colnames
# colnames(human.tf.dat)[2] = 'gsymb'

# ## get mouse genes one way
# #mouse_genes = homologene(human.tf.dat$gsymb,inTax = 9606, outTax = 10090)

# ## get mouse genes second way using a wrapper funciton. COnvenient as it has columns easy to read
# mouse_genes = human2mouse(human.tf.dat$gsymb)

# ## keep selected rows
# mouse_genes = mouse_genes[, c('humanGene','mouseGene')]

# ## convert to upper
# mouse_genes$mouseGene = toupper(mouse_genes$mouseGene)

# ## get DBD info from TF dbd data
# mouse_DBD = inner_join(mouse_genes, human.tf.dat, by = c('humanGene' = 'gsymb'))

# ## split comma separated rows
# mouse_DBD = separate_rows(mouse_DBD, c('DBD'), sep = ';')

# ## strip leading and trailing white spaces
# mouse_DBD$DBD = trimws(mouse_DBD$DBD)

#################################################################################### using CIS-BP ###########################################################
## read DBD file
mouse_DBD = read.delim(paste0('/home/bsharmi6/NA_TF_project/RNA-seq_all_TPM/CIS-BP-2.0.txt'), h = T, stringsAsFactors = F)
## keep relevant columns
mouse_DBD = mouse_DBD[, c('Name', 'Family')]
## convert to upper
mouse_DBD$Name = toupper(mouse_DBD$Name)
## split by comma
mouse_DBD = separate_rows(mouse_DBD, Family, sep = '\\,|\\/')
## change colnames
colnames(mouse_DBD)[2] = 'DBD'

## get common in tpm semi join since we dont want duplicate genes and only TF is needed, not their domains
#tpm.annotated = semi_join(tpm.dat, mouse_DBD, by = c('gene_symbol' = 'mouseGene'))
tpm.annotated = semi_join(tpm.dat, mouse_DBD, by = c('gene_symbol' = 'Name'))

##### mere replicates
##### mere replicates
.replicate.merged <- function(s) {
  Symbol = s[, 1]
  x = s[, 2:ncol(s)]
  
  ## merge replicates
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
  #lab.order = c(8:14, 7, 2, 4:6, 1, 3)  
  #xlab <- xlab[lab.order]; w <- w[,lab.order]
  #xlab = gsub('forebrain_', "", xlab)
  #xlab = gsub('d0', 'P0', xlab)
  #xlab = gsub('22mon', '22mo', xlab)
  ## add gsymb
  w<-cbind(Symbol, w)
  return(list(w,xlab))  
}

## call merge function
res<-.replicate.merged(tpm.annotated)
## get the gexp matrix
gene_exp = res[[1]]
#y_mat = cbind.data.frame(Symbol = y_mat$gsymb, log2(1+y_mat[, 2:ncol(y_mat)]))
xlab = res[[2]]

## remove rows with all 0s
gene_exp = gene_exp[! apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) all(x ==0)),]

## get fold change
#gene_exp = cbind.data.frame(Symbol = gene_exp[, 1], gene_exp[,2:ncol(gene_exp)]/rowMeans(gene_exp[, 2:ncol(gene_exp)], na.rm = T))

######################################################################################### 2 fold change ###############################################################################################################
## set type
#type = c(colnames(gene_exp)[2:ncol(gene_exp)], 'Mixed')

## create a list for cell types
#combined.cluster.celltype = vector('list', length(type))
## set names
#names(combined.cluster.celltype) = type

## get list for each cell tyoe
# for(itype in 1:length(type)-1){
#   ## get row index for max in cell type
#   celltype.max.dat = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) (x[itype] == max(x) & max(x[itype] - x)>=2)),]
#   ## sort
#   celltype.max.dat = celltype.max.dat[order(-celltype.max.dat[, itype+1]),]
# }

## combine all clustering together with a big

######################################################################################### biclustering with BCCC ###############################################################################################################
## get quantile probabilities
#gene_exp_quantile = quantile(unlist(gene_exp[, 2:ncol(gene_exp)]), probs = c(0.15, 0.5, 0.90))

## get index of genes with 0 in more tahn half samples
#gene_exp_zeros = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) sum(x==0)>= length(x)/2),]

## consistent low
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) mean(x, na.rm = T) < gene_exp_quantile[1]),]
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) all(x < gene_exp_quantile[1])),]

## consistent high
#gene_exp_upper_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) mean(x, na.rm = T) > gene_exp_quantile[3]),]
#gene_exp_upper_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) all(x > gene_exp_quantile[3])),]

## write low and high TFs
#write.table(gene_exp_lower_quartile$Symbol, file = paste0(fpath, 'BCCC/forebrain/Low_expressed_TFs.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = F)
#write.table(gene_exp_upper_quartile$Symbol, file = paste0(fpath, 'BCCC/forebrain/High_expressed_TFs.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = F)

## get dynamic genes 
#mixed_genes = setdiff(gene_exp$Symbol, gene_exp_lower_quartile$Symbol)

## middle
#gene_exp_middle_quartile = gene_exp[gene_exp$Symbol %in% mixed_genes, ]

## data matrix must be standardized
# tmp =  t(scale(t(gene_exp[, 2:ncol(gene_exp)])))

# ## perform clustering
# cl <- biclust(as.matrix(tmp), method=BCCC(), delta=0.4,  alpha=2, number=4)

# ## check cl
# cl

# ## check number of clusters
# num_clust = cl@Number

# #### define output list
# moduleLabels_list=vector('list', length(seq(1, num_clust)))
# names(moduleLabels_list) = seq(1, num_clust) 

# ## assign dms to clusters
# for(iclust in 1:length(moduleLabels_list)){
#   ## get row/gene index
#   row_indx = which(cl@RowxNumber[,iclust])
#   ## get column/sample index
#   col_indx = which(cl@NumberxCol[iclust,])
#   moduleLabels_list[[iclust]] = cbind.data.frame(Symbol = gene_exp$Symbol[row_indx], gene_exp[, 2:ncol(gene_exp)][row_indx,col_indx])
# }

# ## check some statistincs
# sapply(moduleLabels_list, function(x) colnames(x))
# sapply(moduleLabels_list, function(x) dim(x))
# setdiff(colnames(gene_exp), unique(unlist(sapply(moduleLabels_list, function(x) colnames(x)))))
# sum(sapply(moduleLabels_list, function(x) dim(x)[1]))

# ## save
# save(cl, file = paste0(fpath, 'Results/BCCC_cluster.RData'))

# ## get gene names in each cluster
# gene_list=vector('list', length(seq(1, num_clust)))
# names(gene_list)=seq(1, num_clust)
# for(iclust in 1:length(gene_list)){
#   ## get row/gene index
#   row_indx = which(cl@RowxNumber[,iclust])
#   ## get genes
#   gene_list[[iclust]]=gene_exp$Symbol[row_indx]
# }

# ## write target dms files 
# sapply(names(gene_list), function (x) write.table(gene_list[[x]], file=paste0(fpath, 'Results/', x, "_genelist.txt"),append = F,quote = F,sep = '\t', row.names = F,col.names = F))

# ## box plots for  different time for each cluster of genes
# for (iclust in 1:length(gene_list)){
#   ##  retrieve genes exp
#   genes.exp.iclust = moduleLabels_list[[iclust]]
#   ## add expression for duplicate gene symbol
#   #genes.exp.iclust = ddply(genes.exp.iclust,"gsymb",numcolwise(mean))
#   ## divide by mean of each gene
#   #genes.exp.iclust = cbind.data.frame(Symbol = genes.exp.iclust[, 1], genes.exp.iclust[,2:ncol(genes.exp.iclust)]/rowMeans(genes.exp.iclust[, 2:ncol(genes.exp.iclust)], na.rm = T))
#   ## change to log
#   genes.exp.iclust = cbind.data.frame(Symbol = genes.exp.iclust[, 1], log(genes.exp.iclust[,2:ncol(genes.exp.iclust)]+1))
#   ## plot
#   png(paste0(fpath, 'Results/BCCC_gene_exp_dev_',names(gene_list)[[iclust]],'.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
#   pheatmap(genes.exp.iclust[, 2:ncol(genes.exp.iclust)], color = cols, cluster_col=F, cluster_rows = F, show_rownames=F, clustering_method="ward.D2", fontsize=24, cellwidth=35)
#   dev.off()
# }

# ############################# results with white color for NA ##################################
# ## combine all clustering together with a big heatmap
# combined.cluster.dat = data.frame(matrix(NA,ncol = ncol(gene_exp)))

# ## set col names
# colnames(combined.cluster.dat) = colnames(gene_exp)

# ## fill up the matrix for first
# for(iclust in c(1)){
# #for(iclust in 1:length(moduleLabels_list)){
#   ## get genes
#   row_indx = which(gene_exp$Symbol %in% gene_list[[iclust]])
#   ## fill up matrix
#   combined.cluster.dat = rbind.data.frame(combined.cluster.dat, gene_exp[row_indx,])
# }

# ## set rownames to null
# rownames(combined.cluster.dat) = NULL

# ## remove first
# combined.cluster.dat = combined.cluster.dat[-1,]

# ## combine all clustering together with a big heatmap
# tmp = data.frame(matrix(NA,ncol = ncol(gene_exp)))

# ## set col names
# colnames(tmp) = colnames(gene_exp)

# ## fill up the matrix for different clusters
# for(iclust in c(2,3,4)){
# #for(iclust in 1:length(moduleLabels_list)){
#   ## get genes
#   row_indx = which(gene_exp$Symbol %in% gene_list[[iclust]])
#   ## fill up matrix
#   tmp = rbind.data.frame(tmp, gene_exp[row_indx,])
# }

# ## remove first
# tmp = tmp[-1,]

# ## sort tmp based on median
# tmp$median_tmp = apply(tmp[, 2:ncol(tmp)], 1, function(x) mean(x, na.rm = T))
# tmp = tmp[order(-tmp$median_tmp),]

# ## set rownames to null
# rownames(tmp) = NULL

# ## remove the column 
# tmp = within(tmp, rm(median_tmp))

# ## combine all clusters
# combined.cluster.dat = rbind.data.frame(combined.cluster.dat, tmp)

# ## add low a upper quartile
# #combined.cluster.dat = rbind.data.frame(combined.cluster.dat, gene_exp_upper_quartile, gene_exp_lower_quartile)


# ## divide by mean of each gene exclusing rows with all 0s
# #combined.cluster.dat = cbind.data.frame(Symbol = combined.cluster.dat[, 1], combined.cluster.dat[,2:ncol(combined.cluster.dat)]/rowMeans(combined.cluster.dat[, 2:ncol(combined.cluster.dat)], na.rm = T))
# combined.cluster.dat = cbind.data.frame(Symbol = combined.cluster.dat[, 1], log(combined.cluster.dat[, 2:ncol(combined.cluster.dat)] +1))

# ## write
# #wite.table(combined.cluster.dat, file = paste0(fpath, 'BCCC/combined.cluster.dat.txt'), append = F,quote = F,sep = '\t', row.names = F,col.names = T)
# ## set NA to 0s
# #na.indx = which(apply(combined.cluster.dat[,2:ncol(combined.cluster.dat)], 1, function(x) all(is.na(x))))

# ## generate combined heatmap
# png(paste0(fpath, 'Results/BCCC_gene_exp_dev_combined_low_mid_high.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
# pheatmap(combined.cluster.dat[, 2:ncol(combined.cluster.dat)], color = cols, cluster_col=F, cluster_rows = F, show_rownames=F, clustering_method="ward.D2", fontsize=24, cellwidth=35)
# dev.off()

## get cluster info
#res <- pheatmap(combined.cluster.dat[, 2:ncol(combined.cluster.dat)], color = cols, cluster_col=F, cluster_rows = T, show_rownames=F, clustering_method="ward.D2", fontsize=24, cellwidth=35)

## get desired groups
#res_groups<- sort(cutree(res$tree_row, k=4))

# ## generate heatmap for lower quartile
# gene_exp_lower_quartile = cbind.data.frame(Symbol = gene_exp_lower_quartile[, 1], log(gene_exp_lower_quartile[,2:ncol(gene_exp_lower_quartile)]+1))
# png(paste0(fpath, 'BCCC/', 'BCCC_gene_exp_dev_combined_low.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
# pheatmap(gene_exp_lower_quartile[, 2:ncol(gene_exp_lower_quartile)], color = cols, cluster_col=F, cluster_rows = F, show_rownames=F, clustering_method="ward.D2", fontsize=24, cellwidth=35)
# dev.off()

# ## generate cheatmap for upper quartile
# gene_exp_upper_quartile = cbind.data.frame(Symbol = gene_exp_upper_quartile[, 1], log(gene_exp_upper_quartile[,2:ncol(gene_exp_upper_quartile)]+1))
# png(paste0(fpath, 'BCCC/', 'BCCC_gene_exp_dev_combined_high.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
# pheatmap(gene_exp_upper_quartile[, 2:ncol(gene_exp_upper_quartile)], color = cols, cluster_col=F, cluster_rows = F, show_rownames=F, clustering_method="ward.D2", fontsize=24, cellwidth=35)
# dev.off()


####################################################################################### annotate the genes for the barplot figure ########################################################################
## get quantile probabilities
gene_exp_quantile = quantile(unlist(gene_exp[, 2:ncol(gene_exp)]), probs = c(0.25, 0.5, 0.7, 0.90))

## consistent low different ways
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) mean(x, na.rm = T) < gene_exp_quantile[1]),]
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) all(x < gene_exp_quantile[1])),]
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) (sum(as.numeric(x)==0)>=ceiling(ncol(gene_exp)*0.5) & quantile(as.numeric(x), probs = 0.9) <= gene_exp_quantile[1])),]
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) (sum(as.numeric(x)==0)>=ceiling(ncol(gene_exp)*0.4) & mean(as.numeric(x)) <= gene_exp_quantile[2] & all(as.numeric(x) <= floor(gene_exp_quantile[3])))),] ## use mean and gene_exp_quantile[2] if including WC sample
gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) (mean(as.numeric(x)) <= gene_exp_quantile[2] & all(as.numeric(x) <= (floor(gene_exp_quantile[3]))))),]
#gene_exp_lower_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) (mean(as.numeric(x)) <= gene_exp_quantile[2] | all(as.numeric(x) <= 1.0))),]

## consistent high
#gene_exp_upper_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) mean(x, na.rm = T) > gene_exp_quantile[3]),]
#gene_exp_upper_quartile = gene_exp[apply(gene_exp[, 2:ncol(gene_exp)], 1, function(x) all(x > gene_exp_quantile[3])),]

## write low and high TFs
#write.table(gene_exp_lower_quartile$Symbol, file = paste0(fpath, 'BCCC/forebrain/Low_expressed_TFs.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = F)
#write.table(gene_exp_upper_quartile$Symbol, file = paste0(fpath, 'BCCC/forebrain/High_expressed_TFs.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = F)

## get dynamic genes 
mixed_genes = setdiff(gene_exp$Symbol, gene_exp_lower_quartile$Symbol)

## middle
gene_exp_middle_quartile = gene_exp[gene_exp$Symbol %in% mixed_genes, ]

## get fold change
gene_exp_middle_quartile = cbind.data.frame(Symbol = gene_exp_middle_quartile[, 1], gene_exp_middle_quartile[,2:ncol(gene_exp_middle_quartile)]/rowMeans(gene_exp_middle_quartile[, 2:ncol(gene_exp_middle_quartile)], na.rm = T))

# set type
type = c(colnames(gene_exp_middle_quartile)[2:ncol(gene_exp_middle_quartile)])

# create a list for cell types
combined.cluster.celltype = vector('list', length(type))
# set names
names(combined.cluster.celltype) = type

# get list for each cell tyoe
for(itype in 1:(length(type))){
  ## get row index for max in cell type
  celltype.max.dat = gene_exp_middle_quartile[apply(gene_exp_middle_quartile[, 2:ncol(gene_exp_middle_quartile)], 1, function(x) (x[itype] == max(x) & min(x[itype] - x[-itype])>=2)),]
  ## sort
  celltype.max.dat = celltype.max.dat[order(-celltype.max.dat[, itype+1]),]
  ## add a column
  celltype.max.dat$Type = rep(type[itype], nrow(celltype.max.dat))
  ## add to list
  combined.cluster.celltype[[itype]] <- celltype.max.dat
}

## cell.type genes
cell.type.dat = do.call('rbind', combined.cluster.celltype)
## set 0 rownames
rownames(cell.type.dat) = NULL

## get the remaiing genes in mixed type
mixed.gene.dat = gene_exp_middle_quartile[gene_exp_middle_quartile$Symbol %in% setdiff(gene_exp_middle_quartile$Symbol, cell.type.dat$Symbol),]

## add type
mixed.gene.dat$Type = rep('Mixed', nrow(mixed.gene.dat))

## sort tmp based on median
mixed.gene.dat$median_tmp = apply(mixed.gene.dat[, 2:(ncol(mixed.gene.dat)-1)], 1, function(x) mean(x, na.rm = T))
mixed.gene.dat = mixed.gene.dat[order(-mixed.gene.dat$median_tmp),]

## remove the column 
mixed.gene.dat = within(mixed.gene.dat, rm(median_tmp))

## convert low expression to log
gene_exp_lower_quartile = cbind.data.frame(Symbol = gene_exp_lower_quartile$Symbol, log2(gene_exp_lower_quartile[, 2:ncol(gene_exp_lower_quartile)] +1))

## sort  based on median
gene_exp_lower_quartile$median_tmp = apply(gene_exp_lower_quartile[, 2:ncol(gene_exp_lower_quartile)], 1, function(x) mean(x, na.rm = T))
gene_exp_lower_quartile = gene_exp_lower_quartile[order(-gene_exp_lower_quartile$median_tmp),]

# ## set rownames to null
rownames(gene_exp_lower_quartile) = NULL

# ## remove the column 
gene_exp_lower_quartile = within(gene_exp_lower_quartile, rm(median_tmp))

## add type for low expressed genes
gene_exp_lower_quartile$Type = rep('Low', nrow(gene_exp_lower_quartile))

## combine all cell type specific, mixed and low
combined.cluster.dat = rbind.data.frame(cell.type.dat, mixed.gene.dat, gene_exp_lower_quartile)

## dchange to log
#combined.cluster.dat = cbind.data.frame(Symbol = combined.cluster.dat[, 1], Type = combined.cluster.dat[,ncol(combined.cluster.dat)], log(combined.cluster.dat[, 2:(ncol(combined.cluster.dat)-1)] +1))

## write
write.table(combined.cluster.dat, file = paste0(fpath, 'Results/combined.cluster.dat.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = T)

# ## generate combined heatmap
png(paste0(fpath, 'Results/BCCC_gene_exp_dev_combined_low_mid_high.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
pheatmap(combined.cluster.dat[, 2:(ncol(combined.cluster.dat)-1)], color = cols, cluster_col=F, cluster_rows = F, show_rownames=F, clustering_method="ward.D2", fontsize=24, cellwidth=35)
dev.off()

## for writing supplementary table
tmp = rbind.data.frame(cell.type.dat, mixed.gene.dat)
## get original tpm and convert to log2
tmp1 = inner_join(tmp, gene_exp, by = 'Symbol')
## log2
tmp1= cbind.data.frame(Symbol = tmp1$Symbol, Type = tmp1$Type, log2(tmp1[, 10:16] +1))
## remove .y
colnames(tmp1) = gsub('.y$', '', colnames(tmp1))
## merge
tmp1 = rbind.data.frame(tmp1, gene_exp_lower_quartile)

## write
write.table(tmp1, file = paste0(fpath, 'Results/Supplementary_Table2.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = T)

############################################################################### plotting ##############################################################

## melt data for ggplot
ggdat <- melt(combined.cluster.dat,id.vars = c("Symbol", "Type"))

# chagne five to 3 groups
ggdat$Type = gsub('High|Low|Mixed', 'Other', ggdat$Type)

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
png(paste0(fpath, 'Results/BCCC_gene_exp_barplot.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
ggplot(ggdat, aes(fill=Type, y=value, x=variable)) + geom_jitter(aes(colour = Type)) +  xlab("") + ylab("Abundance(Log scale)") +
  scale_color_manual(values=cbPalette) +theme_set(theme_classic(base_size = 18)) +theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
dev.off()


################################################################################ get tpm annotated for barplot ###############################
TF.annotated = inner_join(combined.cluster.dat, mouse_DBD, by= c('Symbol' = 'Name'))

## change to 3 groups
TF.annotated$Type = gsub('High|Low|Mixed', 'Other', TF.annotated$Type)

## delete redundant columns
TF.annotated = TF.annotated[, c('Symbol', 'Type', 'DBD')]

## ggdat
ggdat = data.frame(table(TF.annotated$Type, TF.annotated$DBD))

## colnams
colnames(ggdat) = c('Type', 'DBD', 'Count')

## change to 3 groups
#ggdat$Type = gsub('High|Low|Mixed', 'Other', ggdat$Type)

## split by group
ggdat_split = split(ggdat, ggdat$DBD)

## subtract sum by individual
frac_values = lapply(ggdat_split, function(x) (x[,3]/sum(x[,3]))*log(sum(x[,3])+1))

## append to ggdat
ggdat$Count_frac = unlist(frac_values)
write.table(ggdat, file = paste0(fpath, '/Results/ggdat.txt'),append = F,quote = F,sep = '\t', row.names = F,col.names = T)

## convert count to log
#ggdat$Count_log = log2(ggdat$Count+1)
#ggdat = do.call('rbind', ggdat_split)

## define colors
#cbPalette=c("forestgreen","antiquewhite", "antiquewhite1", "pink", 'blue')
cbPalette=c("red","maroon1",  "blue", "chartreuse", "darkorange", "yellowgreen", "lightskyblue2", "khaki")

#ggplot(ggdat, aes(fill=Type, y=Count_frac, x=reorder(DBD, -Count_frac))) + geom_bar(stat="identity", colour = 'black') +  xlab("") + ylab("Number of TFs (log)") + theme_set(theme_classic(base_size = 18)) +
#  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.title=element_blank()) + scale_fill_manual(values=cbPalette) +scale_y_continuous(expand = c(0,0)) 

## plot
png(paste0(fpath, '/Results/BCCC_gene_exp_DBD.png'), type = 'cairo', width = 1800, height = 850); par(mai=c(3,1,1,1))
ggplot(ggdat, aes(fill=Type, y=Count_frac, x=reorder(DBD, -Count_frac))) + geom_bar(stat="identity", colour = 'black') +  xlab("") + ylab("Number of TFs (log)") + theme_set(theme_classic(base_size = 18)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.title=element_blank()) + scale_fill_manual(values=cbPalette) +scale_y_continuous(expand = c(0,0)) 
dev.off()
