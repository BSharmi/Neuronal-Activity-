scRNAseq_analyzer <- function(args=c("/groups/ECBL/project_data/VTCRI_Methylome/scRNAseq.sharmi", "logUMI")) {

  # create an object of Seurat with normalization to 10000 (total counts)
  # reference: https://davetang.org/muse/2017/08/01/getting-started-seurat/
  # install.packages("Seurat", dependencies=TRUE, INSTALL_opts = c('--no-lock'))
  
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(RColorBrewer)
  library(ggplot2)

  ## define output path
  fpath = '/home/bsharmi6/NA_TF_project/scRNAseq_2017_Cell/scRNAseq_analyse/'

  sc.matrix <- c()
  for(dirName in dir(path=args[1], pattern=args[2], full.names=T)) {
    
    if(!dir.exists(dirName)) next
    
    # to load raw data
    sc.rawData <- Read10X(data.dir = dirName, gene.column = 1)
    rownames(sc.rawData) <- paste(unlist(read.table(sprintf("%s/genes.tsv", dirName), h=F)))
    colnames(sc.rawData) <- paste0(colnames(sc.rawData), sub(".*/", "", dirName), sep="")
    sc.rawData <- sc.rawData[order(rownames(sc.rawData)),]
    
    # to merge by time points
    if(is.null(sc.matrix)) sc.matrix <- sc.rawData else sc.matrix <- cbind(sc.matrix, sc.rawData)
  }
  
  # to reset the name to output
  #dirName <- paste(args[1],"/", args[2], sep="")
  dirName <- fpath
  
  # to create Seurat object
  sc <- CreateSeuratObject(counts = sc.matrix, project = "SC",min.cells = 3, min.features = 200)
  
  # mitochondria genes conveniently start with MT
  #mito.genes <- grep(pattern = "^MT-", x = rownames(x = sc@data), value = TRUE, ignore.case=T)
  #percent.mito <- Matrix::colSums(sc@raw.data[mito.genes, ]) / Matrix::colSums(sc@raw.data)
  sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")

  # add some more meta data
  #sc <- AddMetaData(object = sc, metadata = percent.mito, col.name = "percent.mito")
  #VlnPlot(object = sc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  
  ## filter high and low expression genes
  #sc <- FilterCells(object = sc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(500, -Inf), high.thresholds = c(15000, 0.01))
  sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

  # currently there is only log-normalization
  #sc <- NormalizeData(object = sc, normalization.method = "LogNormalize", scale.factor = 10000)
  #hist(colSums(sc@data), breaks = 100, main = "Total expression after normalisation", xlab = "Sum of expression")
  sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # refer to ?FindVariableGenes
  #sc <- FindVariableGenes(object = sc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot=FALSE)
  sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
  
  # build linear model using nUMI and percent.mito
  #sc <- ScaleData(object = sc, vars.to.regress = c("nUMI", "percent.mito"))
  all.genes <- rownames(sc)
  sc <- ScaleData(sc, features = all.genes)

  ## run PCA
  #sc <- RunPCA(object = sc, pc.genes = sc@var.genes, pcs.compute=50, do.print = FALSE)
  sc <- RunPCA(sc, features = VariableFeatures(object = sc),npcs=30)
  # The PrintPCA() function outputs a set of genes that most strongly define a set of principal components.
  #PrintPCA(object = sc, pcs.print = 1:2, genes.print = 5, use.full = FALSE)
  
  bmp(sprintf("%s%s.PCA.png", dirName, 'scRNAanalyze'), type = 'cairo', width=500, height=500); par(mai=rep(1,4))
  DimPlot(object = sc, reduction.use = "pca", label.size=15)
  dev.off()
  
  # the results of the projected PCA can be explored by setting use.full=TRUE in the functions above
  #sc <- ProjectPCA(object = sc, do.print = FALSE)

  # Examine and visualize PCA results a few different ways
  print(sc[["pca"]], dims = 1:5, nfeatures = 5)

  #PCHeatmap(object = sc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
  #PCElbowPlot(object = sc, num.pc=50)
  
  # save.SNN = TRUE saves the SNN so that the clustering algorithm can be rerun using the same graph
  # but with a different resolution value see ?FindClusters
  #sc <- FindClusters(object = sc,  force.recalc=TRUE, genes.use=c("Snap25", "Slc17a6", "Stmn2", "Gad1", "Olig1", "Aqp4", "Cldn5", "Vtn", "Cx3cr1", "Mrc1"), reduction.type = "pca",  dims.use = 1:30, resolution = 0.6, print.output = 0, save.SNN = TRUE)
  #sc <- FindClusters(object = sc,  force.recalc=TRUE, reduction.type = "pca",  dims.use = 1:30, resolution = 0.6, print.output = 0, save.SNN = TRUE)
  sc <- FindNeighbors(sc, dims = 1:30, reduction = "pca", nn.eps = 0.5)
  sc <- FindClusters(sc, resolution = 0.6, n.start = 10)
  
  # visualization
  #sc <- RunTSNE(object = sc, dims.use = 1:30,  do.fast = TRUE)
  sc <- RunUMAP(sc, dims = 1:30, min.dist = 0.75) 
  #sc <- RunTSNE(sc,dims.use = 1:30,max_iter=2000,nthreads = 4)

  # how big are the clusters
  #table(sc@ident)
  
  bmp(sprintf("%s%s.UMAP.png", fpath, 'cluster_by_variablegenes'), type = 'cairo', width=500, height=500); par(mai=rep(1,4))
  #TSNEPlot(object = sc, do.label = TRUE)
  DimPlot(sc, reduction = "umap", pt.size = 0.1) + ggtitle(label = "Visualiation of clusters")
  dev.off()
  
  # find markers for every cluster compared to all remaining cells, report only the positive ones
  #sc.markers <- FindAllMarkers(object = sc, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  sc.markers <- FindAllMarkers(sc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  sc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

  ## top10
  top10 <- sc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
  write.table(top10, file=sprintf("%s%s.top10.markers.txt", dirName, 'scRNAanalyze'), quote=F, sep="\t", row.names=F)
  
  ## save sc
  save(sc, file = paste0(fpath, 'sc.RData'))

  # Check expression of genes of interset.
  #genes_to_check = unique(c("Mog","Mobp" ,"Mbp", "Pecam1", "Cldn5", "Ly6c1", "Aqp4", "Aldh1l1", "Fabp7","Slc1a3", "Snap25","Eno2", "Rbfox3","Pdgfra","Susd5", "Cspg4", "Pdgfrb","Mcam", "Des", "Ascl1", "Dcx", "Dlx2","Gad1","Reln", "Calb1", "Slc17a7", "Neurod6", "Mab21l1", "C1qc", "Cd86", "Ctss", "F13a1", "Aif1",
    #"Snap25", "Slc17a6", "Stmn2", "Gad1", "Olig1", "Aqp4", "Cldn5", "Vtn", "Cx3cr1", "Mrc1", "Bdnf"))
  #genes_to_check = c("Egr1", "Mef2c", "Mef2d", "Mef2a", "Ap4", "Mafa", "Junb", "Tbr1", "Tcf21", "Nkx6.1", "Atoh1", "Neurog2", "Rfx6", "Atf3", "Nur77", 
    #"Fra1", "Tgif2", "Nf1", "Oct4", "Rorgt", "Lhx3", "Sox3")
  genes_to_check = c("Mef2a", "Mef2c", "Mef2d")
  #genes_to_check = c("Lhx1", "Lhx2", "Lhx5", "Lhx6", "Lhx9")
  #genes_to_check = c("Egr1",  "Egr2",  "Egr3",  "Egr4")


  png(sprintf("%s%s%s.VOLIN.png", dirName, 'scRNAanalyze',genes_to_check[1]), type = 'cairo',width=1000, height=750); par(mai=rep(1,4))
  #p <- VlnPlot(object = sc, features = genes_to_check,  point.size.use=0.1); print(p)
  p <- VlnPlot(object = sc, features = genes_to_check); print(p)
  dev.off()
  
  png(sprintf("%s%s%s.Marker.tSNE.png", dirName, 'scRNAanalyze',genes_to_check[1]), type = 'cairo',width=1000, height=750); par(mai=rep(1,4))
  FeaturePlot(object = sc, features = genes_to_check, cols = c("gray", "blue"), reduction = "tsne")
  dev.off()
  
  png(sprintf("%s%s%s.Marker.heatmap.png", dirName, 'scRNAanalyze', genes_to_check[1]), type = 'cairo',width=1500, height=1500); par(mai=rep(1,4))
  p <- DoHeatmap(object = sc, features = top10$gene) +scale_fill_gradientn(colors = c("green", "gray", "red")); print(p)
  dev.off()

  
  png(sprintf("%s%s%s.targetGenes.heatmap.png", dirName, 'scRNAanalyze',genes_to_check[1]), type = 'cairo', width=1500, height=1500); par(mai=rep(1,4))
  #p <- DoHeatmap(object = sc, features = genes_to_check, group.spacing=0.3,
   # col.low = "green", col.mid = "gray", col.high = "red", slim.col.label = TRUE, remove.key = TRUE, cex.col=25, cex.row=25, group.cex=25); print(p)
 # p <- DoHeatmap(object = sc, features = genes_to_check) + scale_fill_gradient(low = "green", mid = "white", high = "red"); print(p)
  DoHeatmap(object = sc, features = genes_to_check) + scale_fill_gradientn(colors = c("green", "gray", "red"))+ theme(text = element_text(size = 20))
  dev.off()

  png(sprintf("%s%s.DotPlot.png", fpath, 'scRNAanalyze'), type = 'cairo', width=500, height=500); par(mai=rep(1,4))
  DotPlot(sc, genes_to_check, plot.legend = T, col.max = 2.5) + coord_flip()
  dev.off()

  # Assign cell type identity to clusters
  # to make transparent of a color
  makeTransparent = function(..., alpha=0.5) {
  
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    alpha = floor(255*alpha)  
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)

    .makeTransparent = function(col, alpha) {
      rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }

    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    return(newColor)
    
  }

  # to calculate the percent of each of cell types 
  percent.ct <- function(cell.types, ident) {
    percent <- c()
    for(ct in cell.types) {
      percent <- rbind(percent, data.frame(celltype=ct, percent=sum(ident %in% ct) / length(ident)))
    }
    percent
  }

  cluster.id <- seq(0, length(table(sc@active.ident))-1)
  cty <- read.table(paste0("/groups/ECBL/project_data/VTCRI_Methylome/scExpression/cell.type.genes (collected).txt"), h=T, sep="\t", stringsAsFactors = F)
  mks <- read.table(paste0('/groups/ECBL/project_data/VTCRI_Methylome/scExpression/PNAS.markers.txt'), h=T, sep="\t", stringsAsFactors= F)  
  markers <- rbind(cty[,c("name", "celltype")], mks[,c("name", "celltype")])
  ## markers from https://github.com/czbiohub/tabula-muris/blob/master/00_data_ingest/02_tissue_analysis_rmd/Brain_Non-Myeloid_facs.Rmd
  tmp = data.frame(name = c("Mobp", "Pecam1", "Cldn5","Susd5", "Cspg4", "Pdgfrb","Mcam", "Des", "Snap25","Eno2", "Rbfox3", "Reln", "Calb1","Slc17a7", "Neurod6", "Mab21l1","Ascl1", "Dcx", "Dlx2"), 
    celltype = c("oligodendrocyte", "endothelial", "endothelial", "oligodendrocyte precursor cells", "oligodendrocyte precursor cells", "Pericyte", "Pericyte", "Pericyte", "pan-neuronal", "pan-neuronal", "pan-neuronal", "Inhibitory neuron", "Inhibitory neuron", "Inhibitory neuron", "Inhibitory neuron", "Inhibitory neuron","NPC","NPC", "NPC"))
  markers <- rbind.data.frame(markers, tmp)
  markers = markers[order(markers$celltype),]
  ## make same names
  markers$celltype = gsub('astrocyte', 'Astrocyte', markers$celltype)
  markers$celltype = gsub('microglia', 'Microglia',markers$celltype)
  markers$celltype = gsub('oligodendrocyte', 'Oligodendrocyte',markers$celltype)
  ## unique
  markers = unique(markers)

  ## get top 10 markers with cell annotation
  top10_list = split(top10, top10$cluster)
  for(itop in 1:length(top10_list)){
    common.genes = intersect(top10_list[[itop]]$gene, markers$name)
    if(length(common.genes)>0){
      top10_list[[itop]]$celltype = unique(markers$celltype[markers$name %in% common.genes])[1]
    }else{
      top10_list[[itop]]$celltype = 'Unkown'
    }
  }

  
  cell.types <- c(unique(markers$celltype), "Unkown")
  n <- length(cell.types)
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:n]
  col.use <- makeTransparent(col_vector, alpha=1)
  names(col.use) <- cell.types
  celltype.id <- unlist(unname(sapply(top10_list, function(x) unique(x$celltype))))
  if(!any(levels(sc@active.ident) %in% cell.types)) sc@active.ident <- plyr::mapvalues(x=sc@active.ident, from=cluster.id, to=celltype.id)
  ## plot
  bmp(sprintf("%s%s.UMAP(mapped).png", dirName, 'scRNAanalyze'), type = 'cairo', width=700, height=650); par(mai=rep(1,4))
  #TSNEPlot(object = sc, do.label = FALSE, colors.use = col.use[celltype.id])
  DimPlot(sc, reduction = "umap", pt.size = 0.1) + ggtitle(label = "Visualiation of clusters")
  dev.off()
  write.table(percent.ct(cell.types, sc@ident), file=sprintf("%s%s.celltype.percent.txt", dirName, 'scRNAanalyze'), quote=F, sep="\t", row.names=F)


  rm(list=ls())
  
}

scRNAseq_analyzer()
