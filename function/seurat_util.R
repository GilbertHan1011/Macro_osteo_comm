#' Function for better to use seurat
#' Author: Gilbert Han
#' Date: 2023.4.6
#' 
#== use harmony to integration-------------------
qcFun <-  function(x){
  x <- PercentageFeatureSet(x, "^mt-", col.name = "percent_mito")
  selected_count <- WhichCells(x, expression =( nCount_RNA > 800 & percent_mito < 20 & nFeature_RNA > 300))
  x <- subset(x, cells = selected_count)
  return(x)
}
runSeurat <- function(x,dim=30){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- FindNeighbors(x, dims = 1:dim)
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:dim)
}
runharmony <- function(x,dim=30){
  require(harmony)
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(x)
  x <- ScaleData(x, features = all.genes)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  x <- RunHarmony(x, "orig.ident",assay.use = "RNA",reduction = "pca")
  x <- FindNeighbors(x, dims = 1:dim,reduction = "harmony")
  x <- FindClusters(x, resolution = 0.5)
  x <- RunUMAP(x, dims = 1:dim,reduction = "harmony")
}
