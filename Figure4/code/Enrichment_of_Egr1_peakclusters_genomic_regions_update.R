rm(list=ls())
library(dplyr)
library(tidyr)

## read table
egr1.peak.ann.dat = read.table('/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/ChIP_seq_analysis/Egr1/ETRM_analysis/egr1/egr1.all.clusters.annotated.txt', h=T)

## split multiple
egr1.peak.ann.dat <- egr1.peak.ann.dat %>% separate_rows(gene.tables, sep = ';')
egr1.peak.ann.dat <- egr1.peak.ann.dat %>% separate_rows(repeat., sep = ';')

## replace loci to category names for CpG
colnames(egr1.peak.ann.dat) <- gsub('CpG', 'CGI', colnames(egr1.peak.ann.dat))
colnames(egr1.peak.ann.dat) <- gsub('CGI.island', 'CGI', colnames(egr1.peak.ann.dat))
egr1.peak.ann.dat$CGI <- sapply(egr1.peak.ann.dat$CGI, function(x) {ifelse (x %in% 'null','null','CGI')})
egr1.peak.ann.dat$CGI.shelves <- sapply(egr1.peak.ann.dat$CGI.shelves, function(x) {ifelse (x %in% 'null','null','CGI.shelves')})
egr1.peak.ann.dat$CGI.shores <- sapply(egr1.peak.ann.dat$CGI.shores, function(x) {ifelse (x %in% 'null','null','CGI.shores')})

## create a column for each value in the categories and fill up with the values
varnames <- c(unique(egr1.peak.ann.dat$gene.tables), unique(egr1.peak.ann.dat$repeat.))
varnames <- varnames[!varnames %in% c('null','Unknown','tRNA','rRNA','scRNA','Other')]
for(ivar in 1:length(varnames)){
  varnames_i = varnames[ivar]
  if(grepl("5pUTR|Intron|Promoter|Exon|3pUTR|Distal_promoter",varnames_i)){
    egr1.peak.ann.dat <- egr1.peak.ann.dat %>% mutate(!!varnames_i := ifelse(gene.tables %in% varnames_i, varnames_i, 'null'))
  }else{
    egr1.peak.ann.dat <- egr1.peak.ann.dat %>% mutate(!!varnames_i := ifelse(repeat. %in% varnames_i, varnames_i, 'null'))
  }  
}

## remove extra columns
egr1.peak.ann.dat <- egr1.peak.ann.dat[, !colnames(egr1.peak.ann.dat) %in% c('gene.tables','repeat.')]

## rename numeric to text
egr1.peak.ann.dat$cluster <- gsub('1', 'one', egr1.peak.ann.dat$cluster)
egr1.peak.ann.dat$cluster <- gsub('2', 'two', egr1.peak.ann.dat$cluster)


################### test of enrichmnet ############
.binom.test <- function(x, N, p, alternative="two.sided", name=NA) {
  test <- c()
  for(i in 1:length(x)) {
    t <- binom.test(x[i], N, p[i], alternative=alternative)
    test <- rbind(test, c(t$null, t$estimate, t$estimate/t$null, t$conf.int[1]/t$null, t$conf.int[2]/t$null, t$p.val))
  }
  colnames(test) <- c("expected", "estimated", "OR", "OR_0.025", "OR_0.975", "p.value")
  rownames(test) <- names(x)
  
  #write.table(data.frame(feature=row.names(test), test), file=sprintf("Enrichment.stat(%s).txt", name), quote=F, sep="\t", row.names=F)
  test
}

.ggplot.bar <- function(pdat, labSize=25, no.names=TRUE, is.single=FALSE, name=NA, x.limits=NA, x.label=NA, fwidth=1000, fheight=500, x.anno=1, y.anno=10, ylim.max=2) {
  
  if(is.na(x.limits) || is.na(x.label)) {
    x.limits <- c("Distal_promoter", "Promoter", "5pUTR", "Exon", "Intron", "3pUTR",
                  "CGI", "CGI.shores", "CGI.shelves", "LINE", "SINE", "LTR", "Low_complexity", "Simple_repeat", "DNA", "Satellite")
    x.label <- c("Distal promoter", "Promoter", "5pUTR", "Exon", "Intron", "3pUTR",
                 "CGI", "CGI.shores", "CGI.shelves", "LINE", "SINE", "LTR", "Low_complexity", "Simple_repeat", "DNA", "Satellite")
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
    png(sprintf('Enrichment(%s, Binomal.test).png', name), type = 'cairo', width=fwidth, height=fheight)
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
    png(sprintf('Enrichment(%s, Binomal.test).png', name), type = 'cairo', width=fwidth, height=fheight)
    par(mfrow=c(1,1), mai=c(1,1,1,1))
    if(no.names) df.dat <- data.frame(feature=row.names(pdat), pdat) else no.names <- pdat
    p <- ggplot(df.dat, aes(x=factor(feature), y=as.double(OR), fill=factor(label))) + 
      geom_bar(stat="identity", position=position_dodge(width=0.95)) +
      labs(title="", x="", y = "Log2 odds ratio") + theme_classic() + scale_fill_discrete(name="Comparison") + 
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

#file <- "DMS.all.annotated.txt"
#x <- read.table(file, h=T)
x = egr1.peak.ann.dat
#idx <- which(grepl("X3.UTR", colnames(x)))
#anno <- list(CD1P23_Math5P23="P23 WT vs. P23 Math5KO", CD1P6_CD1P23="P6 WT vs. P23 WT", CD1P6_Math5P6="P6 WT vs. P6 Math5KO", 
#  Math5P6_Math5P23="P6 Math5KO vs. P23 Math5KO")
anno <- as.list(unique(x$cluster)); names(anno) <-unique(x$cluster)

# single case
for(type in levels(x[,"cluster"])) {
  y <- x[x[,"cluster"] %in% type,5:ncol(x)]
  colnames(y)[1:2] <- c("3'UTR", "5'UTR")
  z <- y!="null"; r <- colSums(z)/nrow(z)
  
  t <- read.table('/home/bsharmi6/Annotation/Distribution.CpG17million.txt', h=T)
  pdat <- .binom.test(colSums(z), nrow(z), t[,2], name=anno[[type]])
  df.dat <- data.frame(feature=rownames(pdat), pdat, label="All")
  .ggplot.bar(df.dat, is.single=T, name=anno[[type]], x.anno=6, y.anno=1.8, ylim.max=2)
} 

# combination
df.dat <- c()
for(type in unique(x$cluster)) {
  y <- x[x[,"cluster"] %in% type,5:ncol(x)]
  #colnames(y)[1:2] <- c("3'UTR", "5'UTR")
  z <- y!="null"; r <- colSums(z)/nrow(z)
  
  t <- read.table('/home/bsharmi6/Annotation/Distribution.CpG17million.txt', h=T)
  pdat <- .binom.test(colSums(z), nrow(z), t[,2], name=type)
  df.dat <- rbind(df.dat, data.frame(feature=rownames(pdat), pdat, label=anno[[type]]))
} 
df.dat[, grep('OR', colnames(df.dat))] <- log2(df.dat[, grep('OR', colnames(df.dat))]+1)
write.table(data.frame(feature=row.names(df.dat), df.dat), file=sprintf("Enrichment.stat(All).txt"), quote=F, sep="\t", row.names=F)

.ggplot.bar(df.dat, labSize=30, no.names=F, is.single=F, name="All", x.anno=6, y.anno=5, ylim.max=6, fwidth=2000, fheight=1300)


