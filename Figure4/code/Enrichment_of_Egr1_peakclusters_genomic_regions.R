## TF motif enrichment
rm(list=ls())

library(pheatmap)
library(plyr)
library(dplyr)
library(tidyverse)
library(stringr)
library(Rtsne)
library(ggplot)

## define fpath 
#fpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMs/Enrichment_motifs_in_neurons/all_motif_locations/DMS/'
fpath = paste0('/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/ChIP_seq_analysis/Egr1/')

## read annotated table. it is annotated with clusters and TF peaks
tmp = read.delim(paste0(fpath, list.files(fpath, pattern = '*annotated.txt')), header = T)

## remove _dms
colnames(tmp) <- gsub('_dms', '', colnames(tmp))

## create a column with cluster annotation
tmp$DMS <- apply(tmp[,4:ncol(tmp)], 1, function(x) paste0(names(which(x!=0)),collapse = ',', sep = ''))

## remove not enriched in any
tmp <- tmp[apply(tmp, 1, function(x) nchar(x['DMS'])>1),]

## rename column
tmp <- rename(tmp, cluster_annotated = DMS)

## split tmp by DMS to get individual rows
tmp <- tmp %>% separate_rows(cluster_annotated, sep = ',')

## rearrange
tmp <- tmp %>% select(chrom, start, end, cluster_annotated,everything())


## set excitatory
#tmp$cluster_annotated = gsub("mDL-2|mDL-1|mDL-3|mL23|mL4|mL6-2|mL6-1|mL5-1|mL5-2|mIn-1", "excitatory", tmp$cluster_annotated)
#tmp$cluster_annotated = gsub("mPv|mSst-1|mNdnf-2|mNdnf-1|mSst-2|mVip", "inhibitory", tmp$cluster_annotated)
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
#colnames(tmp)[1:3] = c('chrom', 'start', 'end')
## remove .motifs
#colnames(tmp) <- gsub('.motifs.HOMER', '', colnames(tmp))
## remove end space
#colnames(tmp) <- gsub('\\.$', '', colnames(tmp))
## capitalize
#colnames(tmp) <- str_to_title(colnames(tmp))
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

## enrichment 
.binom.test <- function(x, N, p, alternative="two.sided", name=NA) {
  test <- c()
  for(i in 1:length(x)) {
    t <- binom.test(x[i], N, p[i], alternative=alternative)
    test <- rbind(test, c(t$null, t$estimate, t$estimate/t$null, t$conf.int[1]/t$null, t$conf.int[2]/t$null, t$p.val))
  }
  colnames(test) <- c("expected", "estimated", "OR", "OR_0.025", "OR_0.975", "p.value")
  rownames(test) <- names(x)
  
  write.table(data.frame(feature=row.names(test), test), file=sprintf("Enrichment.stat(%s).txt", name), quote=F, sep="\t", row.names=F)
  test
}

.ggplot.bar <- function(pdat, labSize=25, no.names=TRUE, is.single=FALSE, name=NA, x.limits=NA, x.label=NA, fwidth=1000, fheight=500, x.anno=1, y.anno=10, ylim.max=2) {
  
  if(is.na(x.limits) || is.na(x.label)) {
    x.limits <- c("Distal_promoter", "Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic",
                  "CGI", "CGI_shore", "CGI_shelf", "LINE", "SINE", "LTR", "Low_complexity", "Simple_repeat", "DNA", "Satellite")
    x.label <- c("Distal promoter", "Promoter", "5'UTR", "Exon", "Intron", "3'UTR", "Intergenic",
                 "CGI", "CGI_shore", "CGI_shelf", "LINE", "SINE", "LTR", "Low_complexity", "Simple_repeat", "DNA", "Satellite")
  }
  
  library(ggplot2)  
  ylim.max <- ylim.max
  labSize <- labSize
  
  label <- function(p.value=as.double(paste(df.dat[,"p.value"]))) {
    labs <- rep("", length(p.value))
    labs[p.value < 0.05] <- "*"
    labs[p.value < 0.01] <- "**"
    labs[p.value < 1e-3] <- "***"
    labs
  }
  
  if(is.single) {
    tiff(sprintf('Enrichment(%s, Binomal.test).tiff', name), width=fwidth, height=fheight)
    par(mfrow=c(1,1), mai=c(1,1,1,1))
    if(no.names) df.dat <- data.frame(feature=row.names(pdat), pdat) else no.names <- pdat
    p <- ggplot(df.dat, aes(x=factor(feature), y=as.double(OR), fill= factor(feature))) + #factor(label, levels=sort(levels(df.dat[,"label"]))))) + 
      geom_bar(stat="identity", position=position_dodge(width=0.95)) +
      labs(title="", x="", y = "Odds ratio") + theme_classic() + scale_fill_discrete(name="Feature") + 
      scale_x_discrete(limits=x.limits, label=x.label) +
      theme(axis.text.y = element_text(size=labSize, color = "black"), axis.text.x = element_text(size=labSize, color = "black", vjust=0.5, hjust=1, angle=90)) + 
      theme(axis.title.y = element_text(size=labSize)) + theme(legend.text = element_text(size = labSize)) + theme(legend.title = element_text(size=labSize)) + 
      theme(legend.position="none", legend.key.size = unit(1, "cm")) + ylim(0,ylim.max) +
      geom_hline(yintercept=1, linetype="dashed", color = "grey30", size=1.2) + 
      geom_text(aes(x=factor(feature), y=as.double(OR), label = label()), position = position_dodge(width=0.95), vjust = -0.1, size = labSize*0.2) +
      annotate(geom="text", x=x.anno, y=y.anno, label="*: p-value<0.05, **: p-value<0.01, ***: p-value<0.001", color="black", size=labSize*0.3)
    print(p)
    dev.off()
  } else {
    tiff(sprintf('Enrichment(%s, Binomal.test).tiff', name), width=fwidth, height=fheight)
    par(mfrow=c(1,1), mai=c(1,1,1,1))
    if(no.names) df.dat <- data.frame(feature=row.names(pdat), pdat) else no.names <- pdat
    p <- ggplot(df.dat, aes(x=factor(feature), y=as.double(OR), fill=factor(label))) + 
      geom_bar(stat="identity", position=position_dodge(width=0.95)) +
      labs(title="", x="", y = "Odds ratio") + theme_classic() + scale_fill_discrete(name="Comparison") + 
      scale_x_discrete(limits=x.limits, label=x.label) +
      theme(axis.text.y = element_text(size=labSize, color = "black"), axis.text.x = element_text(size=labSize, color = "black", vjust=0, hjust=1, angle=90)) + 
      theme(axis.title.y = element_text(size=labSize)) + theme(legend.text = element_text(size = labSize)) + theme(legend.title = element_text(size=labSize)) + 
      theme(legend.position="right", legend.key.size = unit(1, "cm")) + ylim(0,ylim.max) +
      geom_hline(yintercept=1, linetype="dashed", color = "grey30", size=1.2) + 
      geom_text(aes(x=factor(feature), y=as.double(OR), label = label()), position = position_dodge(width=0.95), vjust = -0.1, size = labSize*0.2) +
      annotate(geom="text", x=x.anno, y=y.anno, label="*: p-value<0.05, **: p-value<0.01, ***: p-value<0.001", color="black", size=labSize*0.35)
    print(p)
    dev.off()
  }
  
}

.get.pdat <- function(case, control, name=NA) {
  .getRate <- function(x) {
    idx <- which(grepl("X3.UTR", colnames(x)))
    y <- x[,idx:(ncol(x)-1)]
    colnames(y)[1:2] <- c("3'UTR", "5'UTR")
    
    z <- y!="null"
    list(x=colSums(z), n=nrow(z), rate=colSums(z)/nrow(z))
  }
  
  case.list <- .getRate(read.table(case, h=T))
  control.list <- .getRate(read.table(control, h=T))
  
  .binom.test(case.list$x, case.list$n, control.list$rate, name=name)
}

file <- "DMS.all.annotated.txt"
x <- read.table(file, h=T)
idx <- which(grepl("X3.UTR", colnames(x)))
anno <- list(CD1P23_Math5P23="P23 WT vs. P23 Math5KO", CD1P6_CD1P23="P6 WT vs. P23 WT", CD1P6_Math5P6="P6 WT vs. P6 Math5KO", 
             Math5P6_Math5P23="P6 Math5KO vs. P23 Math5KO")

# single case
for(type in levels(x[,"type"])) {
  y <- x[x[,"type"] %in% type,idx:(ncol(x)-1)]
  colnames(y)[1:2] <- c("3'UTR", "5'UTR")
  z <- y!="null"; r <- colSums(z)/nrow(z)
  
  t <- read.table('Distribution.CpG(17million).txt.gz', h=T)
  pdat <- .binom.test(colSums(z), nrow(z), t[,2], name=anno[[type]])
  df.dat <- data.frame(feature=rownames(pdat), pdat, label="All")
  .ggplot.bar(df.dat, is.single=T, name=anno[[type]], x.anno=6, y.anno=1.8, ylim.max=2)
} 

# combination
df.dat <- c()
for(type in c("CD1P6_CD1P23", "Math5P6_Math5P23", "CD1P6_Math5P6", "CD1P23_Math5P23")) {
  y <- x[x[,"type"] %in% type,idx:(ncol(x)-1)]
  colnames(y)[1:2] <- c("3'UTR", "5'UTR")
  z <- y!="null"; r <- colSums(z)/nrow(z)
  
  t <- read.table('Distribution.CpG(17million).txt.gz', h=T)
  pdat <- .binom.test(colSums(z), nrow(z), t[,2], name=type)
  df.dat <- rbind(df.dat, data.frame(feature=rownames(pdat), pdat, label=anno[[type]]))
} 
write.table(data.frame(feature=row.names(df.dat), df.dat), file=sprintf("Enrichment.stat(All).txt"), quote=F, sep="\t", row.names=F)

.ggplot.bar(df.dat, labSize=30, no.names=F, is.single=F, name="All", x.anno=6, y.anno=1.8, ylim.max=2, fwidth=2000, fheight=1300)

################################################# simple binomial test ############################
dat.pval <- vector('list', length(unique(tmp$cluster_annotated)))
names(dat.pval) <- unique(tmp$cluster_annotated)
dat.est <- vector('list', length(unique(tmp$cluster_annotated)))
names(dat.est) <- unique(tmp$cluster_annotated)
for(itype in 1:length(dat.pval)){
  r <- sum(grepl(names(dat.pval)[itype],tmp$cluster_annotated))
  N <- nrow(tmp)
  binom.res <- binom.test(r,N,0.0625,alternative="greater")
  dat.pval[[itype]] <- binom.res$p.value
  dat.est[[itype]] <- binom.res$estimate
}
