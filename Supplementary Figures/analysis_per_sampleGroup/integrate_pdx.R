library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

setwd("/data/aronow/Kang/SingleCellHBPilot-HB17/integrate/sample_group/")
options(future.globals.maxSize=(120*1024^3)) # set max size as 40G

# 1. read raw data
hb_list = list(); count = 1
dir_filtered = c("hb17-pdx","hb30-pdx")
for (dir in dir_filtered)
{
  if (dir != "../processed/")
  {
    hb.data = Read10X_h5(paste0("../../processed/",dir, "/filtered_feature_bc_matrix.h5"))
    
    hb = CreateSeuratObject(count = hb.data, min.cells = 0, min.features = 0)
    hb[["percent.mt"]] = PercentageFeatureSet(hb, pattern = "^MT-")
    
    hb = hb[, (hb[["nCount_RNA"]]>800) & (hb[["nFeature_RNA"]] > 500) & (hb[["percent.mt"]] < 10)]
    hb = RenameCells(hb, new.names = paste0(dir, "_", colnames(hb)))
    hb$sample = dir
    
    hb_list[count] = hb; count = count + 1
  }
}

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

DefaultAssay(hb.combined) <- "integrated"
# 5. Run the standard workflow for visualization and clustering
hb.combined <- ScaleData(hb.combined, verbose = FALSE)
hb.combined <- RunPCA(hb.combined, npcs = 30, verbose = FALSE)
hb.combined <- RunUMAP(hb.combined, reduction = "pca", dims = 1:30)
hb.combined <- FindNeighbors(hb.combined, reduction = "pca", dims = 1:30)

hb.combined <- FindClusters(hb.combined, resolution = 0.5)
hb.combined <- FindClusters(hb.combined, resolution = 1)
hb.combined <- FindClusters(hb.combined, resolution = 2)

saveRDS(hb.combined, file = "pdx_integrate_reciPCA_perMito10.RDS")
write.table(hb.combined@meta.data, file = "ctype_pdx_perMito10.txt", sep = "\t")
write.table(hb.combined@reductions$umap@cell.embeddings, file = "umap_pdx_perMito10.txt", sep="\t")

# 6. visualization
pdx = readRDS("./pdx_integrate_reciPCA_perMito10.RDS")

# umap of cell class
pdx$Cell.class = "Tumor cells"
DimPlot(pdx, group.by = "Cell.class",cols = c("#66CC00")) + 
        theme(legend.position = c(0.05, 0.87),
              legend.text = element_text(family = "Helvetica", size = 30),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# umap of samples
DimPlot(pdx, group.by = "sample", cols = c("#FF6666","#8E0000")) + 
        theme(legend.position = c(0.05, 0.87),
              legend.text = element_text(family = "Helvetica", size = 30),
              panel.background = element_blank(),
              axis.line = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")
