args <- commandArgs(trailingOnly=TRUE)
######################################## select DMS sites with 10Xs in both TET1KO or TET2KO #########################
## add libraries
library(dplyr)
library(cluster)
library(factoextra)
library(data.table)
####################################################### parameter definition ############################################
### x is the matrix of methylation values on which clustering will be applied
### dms_filt is the set of dms.  It must be a tab ldelimited text file simialr to BED format with three columns for chromosome, start and end columns
#### number of rows of x and dms_filt must match
## data_back is the background file for HOMER. It must be a tab ldelimited text file simialr to BED format with three columns for chromosome, start and end columns
########################################################## read parameter ###############################################
## file path in args[1]
outpath =  args[1] ## path
#### read methylation matrix
x = as.data.frame(fread(args[2]))
## read dms 
dms_filt = as.data.frame(fread(args[3]))
## read background region to prepare for HOMER 
data_back = as.data.frame(fread(args[4]))
################################################################ WGCNA clustering ##########################################
# library(WGCNA)
# options(stringsAsFactors = FALSE)
# allowWGCNAThreads() 
# ## if the number of DMS is >= 100000, increase minModuleSize
# if(dim(x)[1] >= 70000){
# 	net = blockwiseModules(t(x), power = 3, maxBlockSize = 40000, TOMType = "signed", minModuleSize = 10000, reassignThreshold = 1e-6, mergeCutHeight = 0.25, numericLabels = TRUE, networkType = "signed",
# 	replaceMissingAdjacencies = TRUE, pamRespectsDendro = FALSE, loadTOM = FALSE, verbose = 3)
# }else{
# 	net = blockwiseModules(t(x), power = 3, maxBlockSize = 40000, TOMType = "signed", minModuleSize = 5000, reassignThreshold = 1e-6, mergeCutHeight = 0.25, numericLabels = TRUE, networkType = "signed",
# 	replaceMissingAdjacencies = TRUE, pamRespectsDendro = FALSE, loadTOM = FALSE, verbose = 3)
# }
# ## post process
# geneTree = net$dendrograms[[1]]
# moduleLabels = net$colors
# mergedColors = labels2colors(net$colors)
# MEs = net$MEs
# save(net, MEs, moduleLabels, geneTree,file = paste0(outpath,"WGCNA_result.RData"))
# ######## split dms based on cluster 
# ## get modules
# moduleLabels_list=vector('list', length(unique(moduleLabels)))
# names(moduleLabels_list)=unique(moduleLabels)
# for(i in 1:length(moduleLabels_list)){
#   moduleLabels_list[[i]]=dms_filt[moduleLabels %in% names(moduleLabels_list)[i],]
# }
# ## remove unassigned module
# #moduleLabels_list[names(moduleLabels_list)=='0']=NULL
# ## create directories for each cluster
# for(imodule in 1:length(moduleLabels_list)){
# 	## get dirname
# 	dirName = file.path(outpath, paste0('cluster_',names(moduleLabels_list)[imodule]))
# 	## create two directories for missing an non-missing
# 	if(!dir.exists(dirName)) dir.create(dirName)	
# }
# ## write target dms files 
# sapply(names(moduleLabels_list), function (x) write.table(moduleLabels_list[[x]], file=paste0(outpath, paste0('/cluster_',x),"/cluster_" ,x, "_dms.txt"),append = F,quote = F,sep = '\t',
#               row.names = F,col.names = T ))

# ## generate background dms for each target cluster
# rand_back=lapply(moduleLabels_list, function(x) data_back[sample(nrow(data_back), 2*nrow(x)),])
# ## write background files
# sapply(names(rand_back), function (x) write.table(rand_back[[x]], file=paste0(outpath, paste0('/cluster_',x),"/cluster_" ,x, "_background.txt"),append = F,quote = F,sep = '\t',
#               row.names = F,col.names = T ))
####################################################################################### clara clustering ##################
## get optimal cluster number for different subsets
## define output cluster number list
# if(nrow(x)>40000){
# 	niter = 10
# 	nclust = rep(0,niter)
# 	for(isample in 1:niter){
# 		## sample one fourth of data
# 		x_i = sample_n(x, ceiling(dim(x)[1]/3), replace = FALSE, weight = NULL, .env = NULL)
# 		## run optimal cluster test
# 		p=fviz_nbclust(x_i, clara, method = "silhouette", print.summary = TRUE,correct.d=TRUE) +  theme_classic()
# 		## store cluster number
# 		nclust[isample] = which.max(p$data[,2])	
# 	}
# 	## get the number of clusters predicted maximum times
# 	counts <- table(nclust)
# 	final_nclust = as.numeric(names(counts)[which.max(counts)])
# }else{
# 	## run optimal cluster test
# 	p=fviz_nbclust(x, clara, method = "silhouette", print.summary = TRUE,correct.d=TRUE) +  theme_classic()
# 	## get cluster number
# 	final_nclust = which.max(p$data[,2])	
# }

## check cluster number
#p=fviz_nbclust(x, clara, method = "silhouette", print.summary = TRUE,correct.d=TRUE) +  theme_classic()
#final_nclust = which.max(p$data[,2])
# 	## get cluster number
# ## clara
# #clara.res <- clara(x, final_nclust, samples = 5000, pamLike = TRUE, correct.d=TRUE)
# ## clarans
clara.res <- clara(x, 3, samples = 5000, pamLike = FALSE, correct.d=TRUE,rngR = TRUE)
## create CLARA directory
#dir.create(file.path(outpath, paste0('CLARA_results')))	
save(clara.res, file = paste0(outpath, "CLARA_result.RData"))

# ########################################################################### cluster result analysis #########################
# # ## check number of clusters
num_clust = unique(clara.res$cluster)
# #### define output list
moduleLabels_list=vector('list', length(num_clust))
names(moduleLabels_list) = num_clust	
## assign dms to clusters
for(iclust in 1:length(num_clust)){
  clust_indx = which(clara.res$cluster %in% num_clust[iclust])
  moduleLabels_list[[iclust]] = dms_filt[clust_indx,]
}

## create directories for each cluster
for(iclust in 1:length(moduleLabels_list)){
  ## create two directories for missing an non-missing
  dir.create(file.path(outpath, paste0('/cluster_',iclust)))	
}

## write target dms files 
sapply(names(moduleLabels_list), function (x) write.table(moduleLabels_list[[x]], file=paste0(outpath, '/cluster_',x,"/cluster_" ,x, "_dms.txt"),append = F,quote = F,sep = '\t', row.names = F,col.names = T ))

## generate background dms for each target cluster
rand_back=lapply(moduleLabels_list, function(x) data_back[sample(nrow(data_back), 2*nrow(x)),])
## write background files
sapply(names(rand_back), function(x) write.table(rand_back[[x]], file=paste0(outpath, '/cluster_',x ,"/cluster_" ,x, "_background.txt"),append = F,quote = F,sep = '\t', row.names = F,col.names = T ))
