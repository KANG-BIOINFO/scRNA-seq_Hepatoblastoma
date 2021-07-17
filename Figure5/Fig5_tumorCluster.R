# Tumor cell integration (version 2)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

# 0. read tumor cells from all cells
df_tumor = read.table("../sample_group/h5ad/cellxgene_tumor_0217.csv", sep=",", header=1)
df_tumor = df_tumor[df_tumor$Cell.class %in% c("Tumor cell","Neuronal cell"),]
df_pdx = read.table("../sample_group/ctype_pdx_perMito10.txt", sep="\t", header=1, row.names = 1)
tumor_cells = c(as.character(df_tumor$index), rownames(df_pdx))

# 1. separate samples
hb_list = SplitObject(hb_raw_tumor, split.by = "sample")

# 2. normalize and find variable genes in each dataset
hb_list <- lapply(X = hb_list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# 3. select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = hb_list)
hb_list <- lapply(X = hb_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# 4. integrate data
hb.anchors <- FindIntegrationAnchors(object.list = hb_list, anchor.features = features, reduction = "rpca")
hb.combined <- IntegrateData(anchorset = hb.anchors)

# 5. Run the standard workflow for visualization and clustering
DefaultAssay(hb.combined) <- "integrated"
# hb.combined <- ScaleData(hb.combined, vars.to.regress = c("nCount_RNA","nFeature_RNA"), verbose = FALSE)
hb.combined <- ScaleData(hb.combined, verbose = FALSE)
hb.combined <- RunPCA(hb.combined, npcs = 30, verbose = FALSE)
hb.combined <- RunUMAP(hb.combined, reduction = "pca", dims = 1:30)
hb.combined <- FindNeighbors(hb.combined, reduction = "pca", dims = 1:30)

# multiple resolutions
hb.combined <- FindClusters(hb.combined, resolution = 0.2)
hb.combined <- FindClusters(hb.combined, resolution = 0.3)
hb.combined <- FindClusters(hb.combined, resolution = 0.5)
hb.combined <- FindClusters(hb.combined, resolution = 1)
hb.combined <- FindClusters(hb.combined, resolution = 2)

saveRDS(hb.combined, file = "tumor_integrated_v2.RDS")
write.table(hb.combined@meta.data, file = "ctype_tumors.txt", sep="\t")
write.table(hb.combined@reductions$umap@cell.embeddings, file = "umap_tumors.txt", sep="\t")

# 6. visualization
# umap of sample groups of integrated tumor cells (Fig.5)
DimPlot(tumor.integrated, group.by = "integrated_snn_res.0.3",pt.size = 0.2,label.size = 10) +
  theme(legend.position = "none",
        legend.text = element_text(family = "Helvetica", size = 30),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# umap of sample groups of integrated tumor cells (Fig.5)
tumor.integrated$sample_group = factor(tumor.integrated$sample_group, levels = c("tumor","pdx"))
DimPlot(tumor.integrated, group.by = "sample_group",pt.size = 0.2,cols = c("#FFCC99","#FF6666")) +
  theme(legend.position = c(0.7,0.2),
        legend.text = element_text(family = "Helvetica", size = 30),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")


# barplot of sample distribution per cluster (Fig.5)
df_dist = as.data.frame(table(tumor.integrated$sample_group, tumor.integrated$integrated_snn_res.0.3))
colnames(df_dist) = c("sample_group", "cluster", "count")
df_dist$sample_group = factor(df_dist$sample_group, levels = c("tumor", "pdx"))
brks = c(0,0.25,0.5,0.75,1)
ggplot(data = df_dist, aes(x = cluster, y = count, fill = sample_group)) +
  geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, family = "Helvetica",size=30,face="bold"),
        panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),
        axis.line = element_line(),
        axis.title = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm")) +
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks)) +
  scale_fill_manual("sample_group", values = c("tumor" = "#FFCC99", "pdx" = "#FF6666"))

# barplot of sample distribution (Fig.5)
df_dist = as.data.frame(table(tumor.integrated$sample_group, tumor.integrated$integrated_snn_res.0.3))
colnames(df_dist) = c("sample_group", "cluster", "count")
df_dist$sample_group = factor(df_dist$sample_group, levels = c("tumor", "pdx"))
df_dist$cluster = factor(df_dist$cluster, levels = c("Tr-0","Tr-1","Tr-2","Tr-3","Tr-4","Tr-5",
                                                     "Tr-6","Tr-7","Tr-8","Tr-9","Tr-10","Tr-11"))
brks = c(0,0.25,0.5,0.75,1)
# calculate the percentage of cluster in each sample group
df_dist[df_dist$sample_group == "pdx","percent"] = df_dist[df_dist$sample_group == "pdx","count"] / sum(df_dist[df_dist$sample_group == "pdx","count"])
df_dist[df_dist$sample_group == "tumor","percent"] = df_dist[df_dist$sample_group == "tumor","count"] / sum(df_dist[df_dist$sample_group == "tumor","count"])

ggplot(data = df_dist, aes(x = cluster, y = percent, fill = sample_group)) +
  geom_bar(width = 0.9,stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1.0, family = "Helvetica",size=30,face="bold"),
        panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),
        axis.line = element_line(),
        axis.title = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm")) +
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks)) +
  scale_fill_manual("sample_group", values = c("tumor" = "#FFCC99", "pdx" = "#FF6666"))

