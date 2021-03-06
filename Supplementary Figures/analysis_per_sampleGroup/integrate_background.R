library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(plyr)

setwd("/data/aronow/Kang/SingleCellHBPilot-HB17/integrate/sample_group/")
options(future.globals.maxSize=(120*1024^3)) # set max size as 40G

# 1. read raw data
hb_list = list(); count = 1
dir_filtered = c("hb17-background","hb53-background")
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

saveRDS(hb.combined, file = "background_integrate_reciPCA_hb17_hb53_perMito10.RDS")
#saveRDS(hb.combined, file = "background_integrate_reciPCA_hb17_hb53_regressout.RDS")
write.table(hb.combined@meta.data, file = "ctype_background_hb17_hb53_perMito10.txt", sep = "\t")
write.table(hb.combined@reductions$umap@cell.embeddings, file = "umap_background_hb17_hb53_perMito10.txt", sep="\t")

# save as anndata object
SaveH5Seurat(hb.combined, filename = "background_integrate_reciPCA.h5seurat")
Convert("background_integrate_reciPCA.h5seurat", dest = "h5ad")



# 6. visualization
# load cellxgene annotations
hb.combined = readRDS("background_integrate_reciPCA_hb17_hb53_perMito10.RDS")
df_cellxgene = read.table("./h5ad/cellxgene_background_0217.csv", sep = ",", header=1, row.names = 1)
hb.combined$Cell.class = df_cellxgene$Cell.class
hb.combined = hb.combined[,! hb.combined$Cell.class %in% c("Low_quality","Doublet")]
hb.combined$Cell.class_concise = mapvalues(hb.combined$Cell.class,
                                   from = as.character(unique(hb.combined$Cell.class)),
                                   to = c("Hep", "HSC", "Endo", "Kupffer", "B", "Prolif Hep",
                                          "T/NK", "Cho", "Lymphatic"))
write.table(background@meta.data, file = "ctype_background_hb17_hb53_perMito10_v2.txt", sep = "\t")

# 6.0 umap of cell class
DimPlot(hb.combined, group.by = "Cell.class_concise", label = T, label.size = 8) + 
  theme(legend.position = "none",
        legend.text = element_text(family = "Helvetica", size = 45),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")


# 6.1 umap of clusters
DimPlot(hb.combined, group.by = "integrated_snn_res.0.5", label = T, label.size = 8) + 
  theme(legend.position = "none",
        legend.text = element_text(family = "Helvetica", size = 35),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")

# 6.2 umap of samples
hb.combined$sample_concise = mapvalues(hb.combined$sample, from = c("hb17-background","hb53-background"),
                                       to = c("hb17b","hb53b"))
DimPlot(hb.combined, group.by = "sample_concise",cols = c("#009999","#66B2FF")) + 
  theme(legend.position = c(0.05, 0.85),
        legend.text = element_text(family = "Helvetica", size = 35),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA,size=5)) + ggtitle("")


# features
DefaultAssay(hb.combined) = "RNA"
FeaturePlot(hb.combined, features = c("GPC3", "ALB", "FLT1",
                                      "COL3A1","CD68","PTPRC"), cols = c("#D2D2D2","#FF2625"), ncol = 3)

# match with previous annotations
df_ori = read.table("../ctype_original_3_samples_v2.csv", sep = ",", header=1, row.names = 1)
overlap_cells = intersect(rownames(df_ori),colnames(hb.combined))
df_ori = df_ori[overlap_cells,]

hb.combined_overlap = hb.combined[,overlap_cells]
hb.combined_overlap$Cell.class = df_ori$Cell.class

DimPlot(hb.combined_overlap, group.by = "Cell.class", label = T)

# check hb-53b
hb53b = hb.combined[,hb.combined$sample == "hb53-background"]
FeaturePlot(hb53b, features = c("GPC3", "ALB", "FLT1",
                                "COL3A1","CD68","PTPRC"), cols = c("#D2D2D2","#FF2625"), ncol = 3)


