library(CellChat)
library(Seurat)
library(NMF)
library(ggalluvial)
library(ComplexHeatmap)

### preparation ###
# cell chat object for human tumor samples only for investigation of tumor-microenvironment interaction.
HB.integrated = readRDS("../HB_integrate_reciPCA_7samples_harmony_perMito10.RDS")
HB.integrated = HB.integrated[, HB.integrated$sample]
df_cells = read.table("../ctype_integrated_harmony_TumorClustersAdded_cluster11Added.txt",
                      sep = "\t", header = 1, row.names = 1)
overlap_cells = intersect(rownames(df_cells), colnames(HB.integrated))
df_cells = df_cells[overlap_cells, ]; HB.integrated = HB.integrated[,overlap_cells]
# get cells from human tumor samples
HB.integrated@meta.data = df_cells
HB.tumor = HB.integrated[, HB.integrated$sample_group == "tumor"]

# get new cell annotations with tumor clusters and microenvironment cell types
Cell_class_v2 = c()
meta_data = HB.tumor@meta.data
for (i in 1:nrow(meta_data))
{
  if (meta_data[i, 10] == "Tumor cell")
  {Cell_class_v2 = append(Cell_class_v2, as.character(meta_data[i,13]))} else{
    Cell_class_v2 = append(Cell_class_v2, as.character(meta_data[i,10]))}
}
HB.tumor$Cell_class_v2 = Cell_class_v2
HB.tumor = HB.tumor[,HB.tumor$Cell_class_v2 != "10"] # remove cluster 10 due to small cluster size (10 cells)

data.input = GetAssayData(HB.tumor, assay = "RNA", slot = "data")
labels = HB.tumor$Cell_class_v2
meta = data.frame(group = labels, row.names = names(labels))

### create cellchat object ###
cellchat = createCellChat(object = data.input, meta = meta, group.by = "group")

# initiate and preprocessing
cellchat = setIdent(cellchat, ident.use = "group")
CellChatDB = CellChatDB.human # designate the signaling database
cellchat@DB = CellChatDB
cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = projectData(cellchat, PPI.human)

# inference of cell-cell communication network
cellchat = computeCommunProb(cellchat, raw.use = F)
cellchat = computeCommunProbPathway(cellchat)
cellchat = aggregateNet(cellchat)
saveRDS(cellchat,"cellchat_tumor_microenvir_subcluster.RDS")

# aggregated cell-cell communication network
groupSize = as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = "Interaction weights/strength")

# each cell group in the aggregated cell-cell communication network
mat <- cellchat@net$weight
par(mfrow = c(4,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# global pattern
# specific pathway
pathways = c("IGF","EGF","CXCL","CDH","ANGPTL","FGF","BMP")
for (pathway.show in pathways)
{
  vertex.receiver = seq(1,11)
  pdf(paste0("figures_individual/",pathway.show,".pdf"), width = 10, height = 10)
  netVisual_aggregate(cellchat, signaling = pathway.show, vertex.receiver = vertex.receiver, layout = "chord")
  dev.off()
}
pathway.show = c("BMP")
pdf(paste0("figures_individual/",pathway.show,"_contribution.pdf"), width = 5, height = 5)
netAnalysis_contribution(cellchat, signaling = pathway.show)
dev.off()
# netAnalysis_contribution(cellchat, signaling = pathway.show)
# netVisual_aggregate(cellchat, signaling = pathway.show, vertex.receiver = vertex.receiver)

# outgoing
selectK(cellchat, pattern = "outgoing")
nPatterns = 4
par(mfrow = c(1,1), xpd = T)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 5, height = 10, font.size = 3.6)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")

# incoming
selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, width = 5, height = 10, font.size = 3.6)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

# both (grid)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred inte
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

