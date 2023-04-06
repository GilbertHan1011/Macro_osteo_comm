#== read data--------------------
setwd("../11.28_Msx/")
library(colorspace) 
library(RColorBrewer)
library("ggsci")
rm(list=ls())
MsxMerge <- readRDS("2.18_cell_comm/2.18_MsxMerge.Rds")

MsxMerge <- subset(MsxMerge,idents=c("Osteo-lineage", "Perivascular cell", "Macrophages",
                                     "Monocytes","DC", "Neutrophils", "B cells",  "T cells",
                                     "Endothelial cell"))
levels(MsxMerge) <- c("Osteo-lineage", "Perivascular cell", "Macrophages",
                      "Monocytes","DC", "Neutrophils", "B cells",  "T cells",
                      "Endothelial cell")
# assuming your Seurat object is named 'seurat_obj'

# extract the UMAP coordinates of all cells
umap_coords <- as.data.frame(MsxMerge@reductions$umap@cell.embeddings)

# subset based on UMAP coordinates
MsxMerge <- MsxMerge[,umap_coords$UMAP_2 >-14]

# you can also visualize the subsetted cells on UMAP plot
DimPlot(umap_subset, reduction = "umap")
hcl_palettes("diverging", n = 7, plot = TRUE)
DimPlot(MsxMerge,cols = colorspace::qualitative_hcl(n = 10, palette = "Cold"))
DimPlot(MsxMerge,cols = colorRampPalette(brewer.pal(8,'Dark2'))(9))
DimPlot(MsxMerge,cols = pal_jco()(9))
DimPlot(MsxMerge,cols = pal_aaas()(9))
DimPlot(MsxMerge,cols = pal_jama()(9))
DimPlot(MsxMerge,cols = pal_npg()(9))
DimPlot(MsxMerge,cols = pal_lancet()(9))
dir.create("4.5_cell_comm")
ggsave("4.5_cell_comm/4.5_reduction.pdf",width = 10,height = 8)
DimPlot(MsxMerge,cols = pal_npg()(9),label = T)
ggsave("4.5_cell_comm/4.5_reduction.pdf",width = 10,height = 8)
#== cell communication-----------------


library(CellChat)
library(patchwork)
# cellchat <- readRDS("result/2.19_cellchat.Rds")
cellchat <- createCellChat(object = as.matrix(MsxMerge@assays$RNA@counts), meta = as.data.frame(MsxMerge@meta.data) , group.by = "label")
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")) # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
options(future.globals.maxSize= 3891289600)
future::plan("multiprocess", workers = 15) # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("result/2.19_network_visual.pdf")
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength",color.use = pal_npg()(9))
dev.off()
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat2), title.name = rownames(mat)[i])
}
mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat2[4, ] <- mat[4, ]
netVisual_circle(mat2, vertex.weight = groupSize,color.use = pal_npg()(9),
                 weight.scale = T, edge.weight.max = max(mat2), title.name = rownames(mat)[4])


mat3 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
mat3[, 7] <- mat[, 7]
netVisual_circle(mat3, vertex.weight = groupSize, weight.scale = T,color.use = pal_npg()(9),
                 edge.weight.max = max(mat3), title.name = rownames(mat)[7])

#==computer contribute---------------------------------
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
# netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
pathways.show.all <- cellchat@netP$pathways
for (i in pathways.show.all){
  pdf(paste0("result/pathway_show/pathway_",i,".pdf"),width = 10,height = 6)
  netAnalysis_signalingRole_network(cellchat, signaling = i, width = 8, height = 2.5, font.size = 10)
  dev.off()
}
net <- cellchat@net
netp <- cellchat@netP
netProbility <- netp$prob

p <- netAnalysis_signalingRole_network(cellchat, signaling = "GDF", width = 8, height = 2.5, font.size = 10,color.use =pal_npg()(9),color.heatmap = "RdBu")

df.net <- subsetCommunication(cellchat,slot.name = "netP")
test1 <- netProbility[,,1]
test <- apply(netProbility,3,function(x) x[,4])
# calculate column sums
col_sums <- colSums(test)

# subset data frame
test <- test[, col_sums != 0]
pheatmap::pheatmap(test,scale = "row")
test <- test%>%t
library(ComplexHeatmap)
library(circlize)
library(ggsci)
scaled_mat = t(scale(t(test)))
Heatmap(scaled_mat,col=colorRamp2(c(-1, 0, 3),c("#4DBBD5FF", "white", "#E64B35FF")),
        name = "scaled prob",column_title  = "Signaling Pathway Sent by Macrophages")



Heatmap(scaled_mat,col = colorRamp2(c(-1, 0, 2), brewer.pal(3,'YlGn')))

plotGeneExpression(cellchat, signaling = "SPP1")

pathwaySelect <- c("COLLAGEN","FN1","MIF","TGFb","VISFATIN",
                   "GDF","CXCL","SPP1")

outdir <- "result/select_pathway/"
pathways.show="GDF"
for (pathways.show in pathwaySelect){
  pdf(paste0(outdir,"net_",pathways.show,".pdf"),width = 8,height = 6)
  par(mfrow = c(1,1), xpd=TRUE)
  vertex.receiver = seq(1,4) # a numeric vector.
  netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver,color.use =  pal_npg()(9))
  dev.off()
  pdf(paste0(outdir,"net_chord_",pathways.show,".pdf"),width = 8,height = 6)
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
  dev.off()
}
for (pathways.show in pathwaySelect){
  pdf(paste0(outdir,"heatmap_",pathways.show,".pdf"),width = 8,height = 6)
  hm <- netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
  draw(hm)
  dev.off()
}
#> Do heatmap based on a single object
for (pathways.show in pathwaySelect){
  #pdf(paste0(outdir,"contribution_",pathways.show,".pdf"),width = 8,height = 6)
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  ggsave(paste0(outdir,"contribution_",pathways.show,".pdf"))
  #dev.off()
  # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  #pdf(paste0(outdir,"bubble_",pathways.show,".pdf"),width = 8,height = 6)
  netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:7), remove.isolate = FALSE)
  ggsave(paste0(outdir,"bubble_",pathways.show,".pdf"))
  #dev.off()
}
plotGeneExpression(cellchat, signaling = "GDF")
netVisual_chord_gene(cellchat, signaling = c("GDF"),legend.pos.x = 8,color.use = pal_npg()(9))
#> Comparing communications on a single object
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
pdf("result/bubble_macro_heatmap.pdf",width = 10,height = 12)
netVisual_bubble(cellchat, sources.use = 5, targets.use = c(1:10), remove.isolate = FALSE)
dev.off()

pdf("result/bubble_osteo_heatmap.pdf",width = 10,height = 20)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(1:10), remove.isolate = FALSE)
dev.off()
pdf("result/bubble_macro_osteo_heatmap.pdf",width = 10,height = 20)
netVisual_bubble(cellchat, sources.use = c(2,5,6), targets.use = 8, remove.isolate = FALSE)
dev.off()
pdf("result/bubble_osteo_macro_heatmap.pdf",width = 10,height = 20)
netVisual_bubble(cellchat, sources.use = 8, targets.use = c(2,5,6), remove.isolate = FALSE)
dev.off()

pdf("result/chord_osteo_macro_heatmap.pdf",width = 10,height = 8)
netVisual_chord_gene(cellchat, sources.use = 8, targets.use = 5, lab.cex = 0.5,legend.pos.y = 30)
dev.off()
pdf("result/chord_macro_osteo_heatmap.pdf",width = 10,height = 8)
netVisual_chord_gene(cellchat, sources.use = 5, targets.use = 8, lab.cex = 0.5,legend.pos.y = 30)
dev.off()
#> Comparing communications on a single object
#> # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = , targets.use = c(2:7), lab.cex = 0.5,legend.pos.y = 30)



# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathwaySelect)
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2
ggsave("result/signalRole.pdf",width = 10,height = 6)

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",height = 20)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",height = 20)
pdf("result/signalRole_heatmap.pdf",height = 20,width = 10)
ht1 + ht2
dev.off()

library(NMF)
#> Loading required package: pkgmaker
#> Loading required package: registry
#> Loading required package: rngtools
#> Loading required package: cluster
#> NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
#>   To enable shared memory capabilities, try: install.extras('
#> NMF
#> ')
#> 
#> Attaching package: 'NMF'
#> The following objects are masked from 'package:igraph':
#> 
#>     algorithm, compare
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
pdf("result/signal_pattern.pdf",height = 15,width = 15)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns,width = 10,height = 15)
dev.off()


# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
ggsave("result/signal_pattern_river.pdf",width = 10,height = 8)
#> Please make sure you have load `library(ggalluvial)` when running this function
#> 
netAnalysis_dot(cellchat, pattern = "outgoing")
ggsave("result/outgoing_pattern.pdf",width = 10,height = 8)


selectK(cellchat, pattern = "incoming")
nPatterns = 4
pdf("result/signal_pattern_incomming.pdf",height = 15,width = 15)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns,width = 10,height = 15)
dev.off()
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
ggsave("result/incomming_pattern_river.pdf",width = 10,height = 8)
#> Please make sure you have load `library(ggalluvial)` when running this function
netAnalysis_dot(cellchat, pattern = "incoming")
ggsave("result/incoming_pattern.pdf",width = 10,height = 8)


cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
ggsave("result/umap_signal.pdf",width = 8,height = 6)

saveRDS(cellchat,"result/4.6_cellchat.Rds")
