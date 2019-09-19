
library(reshape)
library(ggplot2)
library(stringi)

rm(list=ls())

## read gene exp file
tpm.dat = read.table('/Users/sharmibanerjee/Documents/Fall2018/Methylation/TPM.all.txt', stringsAsFactors = F, sep = '\t', h=T)

## select target genes of interest
#genelist = c('Egr1', 'Junb', 'Tbr1', 'Fosl2')
#genelist = c('Mef2a', 'Mef2b', 'Mef2c', 'Mef2d')
#genelist = c('Lhx3', 'Tfap4', 'Mafa', 'Tcf21', 'Atf3', 'Tbr1','Nr4a1','Pou5f1','Egr1','Rfx2','Junb','Neurog2','Atoh1',
#             'Rorgt','Nf1','Fosl1','Mef2b')
#genelist = c('Tgif2', 'Nr4a1', 'Batf', 'Bach2')
#genelist = c('Sox3', 'Sox10', 'Sox2', 'Sox6', 'Sox4', 'Sox9', 'Sox15', 'Sox17')
## methyl minus
#genelist = c('Atf3', 'Atf1', 'Atf2', 'Junb')
#genelist = c('Zic3', 'Zic1', 'Scl', 'Tcf12')
#genelist = c('Nf1', 'Tlx', 'Tgif2', 'Tbx5', 'Atoh1', 'Rfx6', 'Ascl1', 'Olig2', 'Ar-half', 'Neurog2', 'Ptf1a', 'Neurod1')
#genelist = c('Zeb1', 'E2a', 'Slug', 'Ptf1a', 'Ascl1', 'Heb')
#genelist =c('Egr1','Junb','Mef2c','Mef2d','Mef2c','Nr4a1')
#genelist =c('Tgif2','Sox3','Atf3')
genelist =c('Egr1', 'Egr2', 'Egr3', 'Egr4', 'Mef2c')

## retrieve rows from gene exp
genes.exp = tpm.dat[grepl(paste(genelist, collapse="|"), x = tpm.dat$Symbol, ignore.case = F),]
## another way to retrieve genes exp
genes.exp= tpm.dat[ tpm.dat$Symbol %in% genelist,]

## change to caps
genes.exp$Symbol = gsub('Nf1', 'NF1', genes.exp$Symbol, ignore.case = T)
genes.exp$Symbol = gsub('Nr4a1', 'Nur77', genes.exp$Symbol, ignore.case = T)
genes.exp$Symbol = gsub('Pou5f1', 'Oct', genes.exp$Symbol, ignore.case = T)
genes.exp$Symbol = gsub('Fosl1', 'Fra1', genes.exp$Symbol, ignore.case = T)
genes.exp$Symbol = gsub('Tfap4', 'Ap4', genes.exp$Symbol, ignore.case = T)
## remove some columns
genes.exp = genes.exp[, !grepl(paste(c('neurons', 'hind', 'mid', 'Fetal'), collapse="|"), colnames(genes.exp),ignore.case = T)]
rownames(genes.exp)=genes.exp$Symbol

##### mere replicates
.replicate.merged <- function(x) {
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
  return(list(w,xlab))
  
}


## call merge function
res<-.replicate.merged(genes.exp[,2:ncol(genes.exp)])
tmp = res[[1]]
tmp = log2(1+tmp)
xlab = res[[2]]
## plot
tmp$gsymb = rownames(tmp)
ggdat= melt(tmp,id.vars = 'gsymb')
ggdat = ggdat[order(ggdat$gsymb),]
colnames(ggdat) = c('Genes', 'Time', 'TPM')
ggdat = ggdat[, c("Time", "TPM", "Genes")]
## plot
#ggplot() + geom_line(aes(y = value, x = variable, colour = gsymb),data = ggdat, stat="identity")
## individual plot
#coltype = c('olivedrab1','red', 'royalblue1', 'darkgreen', 'darkmagenta', 'aquamarine','palevioletred1','sienna3')
## TET1KO NF1 
#coltype = c('mediumseagreen','yellow', 'orange3', 'purple1', 'orange3', 'orange3', 'orange3', 'orange3','moccasin', 'orange3')
## TET1KO Mef2 
#coltype = c('royalblue1','red', 'darkgreen', 'green', 'darkorchid2')

## TET1KO Egr 
#coltype = c('royalblue1','red', 'darkgreen')

## TET1KO Sox 
#coltype = rep('red',8)

## TET1KO Lhx 
#coltype = rep('moccasin',6)

## general 
coltype = c('olivedrab1','red', 'royalblue1', 'darkgreen', 'darkmagenta', 'aquamarine','palevioletred1',
            'sienna3', 'black', 'goldenrod2', 'darkslategray1', 'darkorchid1','darksalmon','darkseagreen2',
            'deeppink1','gray48')

par(mai=c(1,1,1,1))
plot(c(1:14), tmp[1,1:(ncol(tmp)-1)], type = 'o', pch = 16, xaxt = "n", xlab="", ylab="log2(1+TPM)", col=coltype[1], lwd=2, ylim = c(min(ggdat$TPM),max(ggdat$TPM)), cex.lab=1.2, cex.axis=1.2)
#text(locator(), labels = rownames(tmp)[1],cex=1.4)
axis(1, at=1:14, labels=xlab, las =2, cex.axis = 1.2)
for (iplot in 2:nrow(tmp)){
  lines(c(1:14), tmp[iplot,1:(ncol(tmp)-1)], type = 'o',pch = 16, col = coltype[iplot], lwd=2)
  #text(locator(), labels = rownames(tmp)[iplot],cex=1.4)
}

#text(locator(), labels = c(unique(ggdat$Genes)))
# Add a legend
#legend("topleft",legend=rownames(tmp),col=coltype, lty=1, cex=0.8, lwd=2, xpd = TRUE, horiz = F, inset = c(0, 0))
#legend("topleft",legend=rownames(tmp),col=coltype, lty=1, cex=1.2, lwd=2,bty = "n",y.intersp=1.2)
#legend("topright",legend=rownames(tmp),col=coltype, lty=1, cex=1.2, lwd=2,bty = "n",y.intersp=1.2)
#legend(11,8.5,legend=rownames(tmp),col=coltype, lty=1, cex=1.2, lwd=2,bty = "n",y.intersp=1.2)
##for mef2
legend(9.5,5.5,legend=rownames(tmp),col=coltype, lty=1, pch = 16, cex=1.2, lwd=2,bty = "n",y.intersp=0.8)

## plot to a file
#png('excitatory.png', type="cairo", width = 600, height = 500)
par(mai=c(1,1,0.5,0.5))
plot(c(1:14), tmp[1,1:(ncol(tmp)-1)], type = 'o', pch = 16, xaxt = "n", xlab="", ylab="log2(1+TPM)", col=coltype[1], lwd=2, ylim = c(min(ggdat$TPM),max(ggdat$TPM)), cex.lab=1.8, cex.axis=1.8)
#text(locator(), labels = rownames(tmp)[1],cex=1.4)
axis(1, at=1:14, labels=xlab, las =2, cex.axis = 1.8)
for (iplot in 2:nrow(tmp)){
  lines(c(1:14), tmp[iplot,1:(ncol(tmp)-1)], type = 'o',pch = 16, col = coltype[iplot], lwd=2)
  #text(locator(), labels = rownames(tmp)[iplot],cex=1.4)
}
## egr1
#legend("bottomright",legend=rownames(tmp),col=coltype, lty=1, pch = 16, cex=1.8, lwd=2,bty = "n",y.intersp=0.8)
## mef2
legend(9.3,2.5,legend=rownames(tmp),col=coltype, lty=1, pch = 16, cex=1.5, lwd=2,bty = "n",y.intersp=0.8)
## bach2
#legend(10,5.5,legend=rownames(tmp),col=coltype, lty=1, pch = 16, cex=1.8, lwd=2,bty = "n",y.intersp=0.8)
#dev.off()
################################################################## methylation plot ###########################################
rm(list=ls())
load('/Users/sharmibanerjee/Documents/Fall2018/Methylation/TEt1ko_TET2ko/TPM_analysis/lhx_TRM_individual_avg_meth.RData')

## merge result sfor ggplot
ggdat = do.call(rbind, final_meth_list)
#ggplot(ggdat, aes(Distance, Frequency, colour=Category)) + geom_line() + theme_classic()
## individual plot
#coltype = c('olivedrab1','red', 'royalblue1', 'darkgreen', 'darkmagenta', 'aquamarine','palevioletred1','sienna3')
## atf
#coltype = c('olivedrab1','royalblue1', 'red', 'darkgreen', 'red', 'darkmagenta','palevioletred1','sienna3')
## lhx
coltype = c('aquamarine','royalblue1', 'olivedrab1', 'darkgreen', 'red', 'darkmagenta','palevioletred1','sienna3')
names(final_meth_list) = gsub(pattern = "nkx6.", "Nkx6-1", names(final_meth_list))
## mef2
#coltype = c('royalblue1','darkgreen', 'red', 'olivedrab1', 'red', 'darkmagenta','palevioletred1','sienna3')
## egr1
#coltype = c('red','olivedrab1', 'royalblue1', 'darkmagenta', 'olivedrab1', 'darkmagenta','aquamarine','sienna3')
## sox
#coltype = c('royalblue1','darkgreen', 'red', 'palevioletred1', 'darkmagenta', 'sienna3','aquamarine','olivedrab1')

plot(c(1:14), final_meth_list[[1]], type = 'o', pch = 16, xaxt = "n", xlab="", ylab="Methylation level", col=coltype[1], lwd=2, ylim = c(min(unlist(final_meth_list)),max(unlist(final_meth_list))), 
     cex.lab=1.2, cex.axis=1.2)
axis(1, at=1:14, labels=names(final_meth_list[[1]]), las =2, cex.axis = 1.2)
for (iplot in 2:length(final_meth_list)){
  lines(c(1:14), final_meth_list[[iplot]], type = 'o',pch = 16, col = coltype[iplot], lwd=2)
}
# Add a legend
#legend("topleft", legend=stri_trans_totitle(names(final_meth_list)),col=coltype, lty=1, cex=1.2, lwd=2,bty = "n",y.intersp=1.2)
legend(9,0.81, legend=stri_trans_totitle(names(final_meth_list)),col=coltype, lty=1, cex=1.2, lwd=2,bty = "n",y.intersp=0.8)



