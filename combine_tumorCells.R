library(Seurat)
library(dplyr)

setwd("/data/aronow/Kang/SingleCellHBPilot-HB17/results/cellxgene/")

# read
tumor <- readRDS("./tumor_v2.RDS")
pdx <- readRDS("./PDX_v2.RDS")

# merge data
tumor$tissue <- "HB T"
pdx$tissue <- "PDX"
tumor$seurat_clusters_individual <- paste("HB T-",tumor$seurat_clusters,sep="")
pdx$seurat_clusters_individual <- paste("PDX-",pdx$seurat_clusters,sep="")
tumor_subset <- tumor[,tumor$cell_type=="Tumor cells"]

reference.list <- c(tumor_subset,pdx)
anchors <- FindIntegrationAnchors(object.list = reference.list,dim=1:30)

merged.tumors <- IntegrateData(anchorset = anchors,dim=1:30)
merged.tumors <- ScaleData(merged.tumors, verbose = F)
merged.tumors <- RunPCA(merged.tumors, npcs = 30, verbose = F)
merged.tumors <- RunUMAP(merged.tumors, reduction = "pca", dims=1:30)

# clustering
merged.tumors <- FindNeighbors(merged.tumors,dim=1:30)
merged.tumors <- FindClusters(merged.tumors,resolution = 0.5)

DimPlot(merged.tumors)

saveRDS(merged.tumors, file = "merged.tumors_v1.RDS")


