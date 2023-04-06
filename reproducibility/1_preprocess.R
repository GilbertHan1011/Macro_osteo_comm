#'reproducibility of preprocess with Seurat
#'Author: Gilbert Han
#'Date: 2023.4.6

#== environment--------
setwd("../11.28_Msx/")

library(Seurat)
library(SingleR)

source("../function/seurat_util.R")


#== read data-------------------
file <- list.dirs("data/",full.names = T,recursive = F)


count10x <- lapply(file,Read10X)
fileName <- file%>%
  gsub("data//","",.)
MsxList <- map2(count10x,fileName,function(x,y) CreateSeuratObject(x,min.cells = 3, min.features = 500,project = y))
lapply(MsxList, dim)

#== quality control----------------
MsxList <- lapply(MsxList,qcFun)
lapply(MsxList, dim)

MsxMerge <- merge(MsxList[[1]],MsxList[2:6])

#== merge------------------
MsxMerge <- runharmony(MsxMerge)

#== annotation---------------------------
DimPlot(MsxMerge)
FeaturePlot(MsxMerge,c("Ptprc","Sp7","Aspn","Adgre1"))
FeaturePlot(MsxMerge,c("Col1a1","Cdh5","Acta2","Cd14","Pdgfrb","Ptprc"),label = T)
Idents(MsxMerge) <- MsxMerge$seurat_clusters
# new.id <- c("Immune-lineage", "Immune-lineage", "Immune-lineage", "Immune-lineage", 
#             "Osteo-lineage", "Osteo-lineage", "Immune-lineage", 
#             "Immune-lineage", "Perivascular cell", "Endothelial cell", "Immune-lineage", "Immune-lineage", 
#             "Immune-lineage", "Osteo-lineage", "Immune-lineage", "Immune-lineage", "Immune-lineage")
new.id <- c("Immune-lineage", "Osteo-lineage", "Immune-lineage", "Immune-lineage", 
            "Immune-lineage", "Immune-lineage", "Immune-lineage", "Immune-lineage", 
            "Perivascular cell", "Endothelial cell", "Immune-lineage", "Immune-lineage", 
            "Osteo-lineage", "Immune-lineage", "Osteo-lineage", "Immune-lineage", "Immune-lineage", "Immune-lineage")
names(new.id) <- levels(MsxMerge)
MsxMerge <- RenameIdents(MsxMerge,new.id)


FeaturePlot(MsxMerge,c("Adgre1","Mrc1","Fcer2a",'Ccl2',"Chil3"))
library(celldex)

ref <- BlueprintEncodeData()
ref <- ImmGenData()
msxSce <- as.SingleCellExperiment(MsxMerge)
pred.hesc <- SingleR(test = msxSce, ref = ref, assay.type.test=1,
                     labels = ref$label.main)
table(pred.hesc$labels)
MsxMerge$pred <- pred.hesc$labels
DimPlot(MsxMerge,group.by = "pred")
plotScoreHeatmap(pred.hesc)
plotDeltaDistribution(pred.hesc, ncol = 3)

ImmuneCell <- subset(MsxMerge,idents=c("Immune-lineage"))
immuneSce <- as.SingleCellExperiment(ImmuneCell)
pred.hesc <- SingleR(test = immuneSce, ref = ref, assay.type.test=1,
                     labels = ref$label.main)
# table(pred.hesc$labels)
# MsxMerge$pred <- pred.hesc$labels
# DimPlot(MsxMerge,group.by = "pred")
ImmuneCell$pred <- pred.hesc$labels
ImmuneCell$label_fine <- pred.hesc$pruned.labels
DimPlot(ImmuneCell,group.by = "pred")
Idents(ImmuneCell) <- ImmuneCell$seurat_clusters
DimPlot(ImmuneCell)
ImmuneCell <- RunTSNE(ImmuneCell)
ImmuneCell <- FindNeighbors(ImmuneCell,reduction = "harmony",dims = 1:20)
ImmuneCell <- Seurat::FindClusters(ImmuneCell)
new.id <- c("Macrophages", "Macrophages", "Monocytes", "Macrophages", 
            "T cells", "T cells", "T cells", "DC", 
            "Neutrophils", "Macrophages", "B cells", "Macrophages", 
            "NKT", "DC", "Macrophages", "Macrophages", "B cells", "Fibroblasts")
new.id <- c("Macrophages", "T cells", "Macrophages", "Macrophages", "Monocytes", "T cells", "B cells", "Neutrophils",
            "Macrophages", "T cells", "DC", "Macrophages", 
            "B cells", "B cells", "DC", "Macrophages", "Macrophages", "Macrophages", "Fibroblasts", "T cells")
names(new.id) <- levels(ImmuneCell)
ImmuneCell <- RenameIdents(ImmuneCell,new.id)
ImmuneCell$label <- Idents(ImmuneCell)
DimPlot(ImmuneCell)
MsxMerge$label <- as.character(Idents(MsxMerge))
MsxMerge_bk <- MsxMerge
MsxMerge[,colnames(ImmuneCell)] <- ImmuneCell
label <- MsxMerge$label
label[colnames(ImmuneCell)] <- as.character(ImmuneCell$label)
MsxMerge$label <- label
Idents(MsxMerge) <- MsxMerge$label
DimPlot(MsxMerge)
saveRDS(MsxMerge,"2.18_MsxMerge.Rds")
DimPlot(MsxMerge,split.by = "orig.ident")


Idents(MsxMerge) <- MsxMerge$orig.ident
new.id <- c("W1", "W1", "W2", "W2", "W2", "W2")
names(new.id) <- levels(MsxMerge)
MsxMerge <- RenameIdents(MsxMerge,new.id)
MsxMerge$time <- Idents(MsxMerge)
Idents(MsxMerge) <- MsxMerge$orig.ident
new.id <-c("C", "N", "C", "C", "N", "N")
names(new.id) <- levels(MsxMerge)
MsxMerge <- RenameIdents(MsxMerge,new.id)
MsxMerge$condition <- Idents(MsxMerge)
Idents(MsxMerge) <- MsxMerge$label
DimPlot(MsxMerge,label = T)
ggsave("result/overall_umap.pdf",width = 8,height = 6)


DimPlot(MsxMerge,label = T,split.by  = "time")
ggsave("result/time_umap.pdf",width = 10,height = 6)
DimPlot(MsxMerge,label = T,split.by  = "condition")
ggsave("result/time_umap.pdf",width = 10,height = 6)


