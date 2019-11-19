library(plyr)
library(dplyr)
library(readxl)
library(pheatmap)
library(RColorBrewer)
library(tidyr)
library(purrr)
library(foreach)
library(doParallel)
library(fst)
library(data.table)

## parallel processing
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
#cl <- makeCluster(20)
registerDoParallel(cl)

## define path
fpath ='/home/bsharmi6/NA_TF_project/scMethylome/'

############################################## read the methylome files for dev and cell type
## read oligo and astro
cell.type.CpG = fread('gunzip -c /groups/ECBL/mouse_brain/mouse_brain/Celltype.ML.ge5X.txt.gz', header = T, stringsAsFactors = F)
## read development data
dev.CpG = fread('gunzip -c /groups/ECBL/mouse_brain/mouse_brain/Dev.ML.ge5X.txt.gz', header = T, stringsAsFactors = F)
## full join two tables
tmp1 = inner_join(dev.CpG, cell.type.CpG, by=c('chrom', 'position'))
## remove high memory tables
rm(cell.type.CpG, dev.CpG)
## add encocde data
encode.CpG = fread('gunzip -c /groups/ECBL/mouse_brain/mouse_brain/ENCODE.ML.ge5X.txt.gz', header = T, stringsAsFactors = F)
## full join two tables
bulk.CpG = inner_join(encode.CpG, tmp1, by=c('chrom', 'position'))
## rename cols for later matching
colnames(bulk.CpG)[1:2] <- c('chr', 'pos')
## remove high memory tables
rm(encode.CpG, tmp1)

## number of columns of output matrix
colnames_output = colnames(bulk.CpG)[3:ncol(bulk.CpG)]
## compress
path <- paste0(fpath, "motif_multi.fst")
write_fst(bulk.CpG, path)
Cpg <- fst(path)
## delete large df
rm(bulk.CpG)

get_methylation_in_TF_peaks <- function(dat1_ipeak, Cpg_f){
	##check for bulk CpG overlap
	CpG.overlap.indx = which(Cpg_f$chr %in% dat1_ipeak$chrom & Cpg_f$pos >= dat1_ipeak$start & Cpg_f$pos <= dat1_ipeak$end)
	if(length(CpG.overlap.indx)>1){
		CpG.overlap = colMeans(Cpg_f[CpG.overlap.indx, 3:ncol(Cpg_f)],na.rm = T)
    }else if(length(CpG.overlap.indx) == 0){
    	CpG.overlap = rep(NA, length(3:ncol(Cpg_f)))
    }else{
    	CpG.overlap = Cpg_f[CpG.overlap.indx, 3:ncol(Cpg_f)]
 	  }   
   ## add to output
   return(CpG.overlap)
}

############################### getting peaks for egr1 motif sites in 16 dmrs #######################
## define TF
TF = 'Egr'
## read file
tf.dat <- read.table(paste0('/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/Enrichment_motifs_in_neurons/ETRM_analysis/', TF, '/', TF, '.trm.txt'), h=T)

##### parallel for loop 
res <- foreach (ipeak = 1 : nrow(tf.dat),.combine = 'list',.multicombine = TRUE,.maxcombine = nrow(tf.dat),.verbose=TRUE, .packages=c("doParallel", "fst")) %dopar% {
   get_methylation_in_TF_peaks(tf.dat[ipeak,], Cpg)
}
##save
save(res, file = paste0(fpath, '/ETRMS_dev_cell_sc/Enrichment_motifs_in_neurons/ETRM_analysis/',TF, '/', TF,'_res_fst.RData'))

## unlist res and append to output dat
res_unlist = data.frame(matrix(unlist(res), nrow=length(res), byrow=T))
out.dat <- cbind.data.frame(tf.dat, res_unlist)
colnames(out.dat)[4:ncol(out.dat)] <- colnames_output
## save
save(out.dat, file = paste0(fpath, '/ETRMS_dev_cell_sc/Enrichment_motifs_in_neurons/ETRM_analysis/',TF,'/',TF,'_peaks_methylation.RData'))

## delete fst object
unlink(path)

################################### create methylaiton matrix files
#rm(list=ls())
#outpath = '/home/bsharmi6/NA_TF_project/scMethylome/ETRMS_dev_cell_sc/'
#filelist = list.files(outpath, pattern = '^[m*]', ignore.case = T)
#for (icell in 1:length(filelist)){
#	meth.mat <- get(load(paste0(outpath, filelist[icell], '/', filelist[icell],'_peaks_methylation.RData')))
#	meth.mat <- meth.mat[, 4:ncol(meth.mat)]
#	write.table(meth.mat, file = paste0(outpath, filelist[icell], '/', filelist[icell], '_methylation.matrix.txt'), append = F,quote = F,sep = '\t', row.names = F,col.names = T )
#}



############ temmporary code #######
# rm(list=ls())

# library(readxl)
# ## define path
# fpath ='/home/bsharmi6/NA_TF_project/scMethylome/'

# for (icell in 1:16){
# 	## read cluster
# 	scmethylome.df = read_excel(paste0(fpath,'Science-table-S5_DMR.xlsx'), sheet = icell, skip = 2)
# 	## get name
# 	scname = as.character(unique(scmethylome.df[,4]))
# 	## get dirname
# 	dirName = paste0(fpath, 'ETRMs/',scname)
# 	## create dir if doesnt exist
# 	if(!dir.exists(dirName)) dir.create(dirName)
# 	## move the DMS, background and methylation matrix file
# 	file.copy(paste0(fpath, 'ETRMs/old_analysis/', scname, "/", scname, "_background.txt"), paste0(fpath, 'ETRMs/', scname))
# 	file.copy(paste0(fpath, 'ETRMs/old_analysis/', scname, "/", scname, "_dms.txt"), paste0(fpath, 'ETRMs/', scname))
# 	file.copy(paste0(fpath, 'ETRMs/old_analysis/', scname, "/", scname, "_methylation.matrix.txt"), paste0(fpath, 'ETRMs/', scname))
# 	file.copy(paste0(fpath, 'ETRMs/old_analysis/', scname, "/", scname, "_peaks_methylation.RData"), paste0(fpath, 'ETRMs/', scname))
# }
