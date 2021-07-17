library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(ggpubr)

# load data
hb.integrated = readRDS("./HB_integrate_reciPCA_7samples_harmony_perMito10.RDS") # all cells
hb.3groups = readRDS("HB_integrate_reciPCA_7samples_harmony_perMito10_3types.RDS") # all hepatocytes and tumor cells
tumor.integrated = readRDS("./tumor integration/tumor_integrated_v3.RDS") # integrated tumor cells

# barplot of cell types in 2 sample_group (Fig.S1)
df_dist = as.data.frame(table(hb.integrated$sample_group, hb.integrated$Cell.class))
colnames(df_dist) = c("sample_group","Cell_class","count")
df_dist =  df_dist[df_dist$count != 0,]
df_dist = df_dist[df_dist$sample_group != "pdx",]
ggplot(data = df_dist, aes(x = sample_group, y = count, fill = Cell_class)) + geom_bar(stat = "identity", position = "fill") +
  theme(axis.text.x = element_text(family = "Helvetica",size=20,face="bold"),
        panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),
        axis.line = element_line(),
        axis.title = element_blank(),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12),
        #legend.position = "none",
        plot.margin = unit(c(1,1,1,1),"cm")) +
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks))

# Violin plot of sequencing depth (Fig.S3b)
VlnPlot(tumor.integrated, group.by = "integrated_snn_res.0.3", features = c("nCount_RNA"),
        pt.size = 0, y.max = 20000) +
  theme(legend.position = "none",
        axis.text = element_text(family = "Helvetica", size = 25),
        axis.title = element_blank())

# tumor cluster enrichment (Fig.S3d)
enrich_file = read.table("./tumor integration/tumor_cluster_enrichment.txt", sep = "\t", header = 1)
tr0 = enrich_file[enrich_file$Cluster == "Tr0",]
barplot(tr0$Score, horiz = T, col = "#FF6666")

tr2 = enrich_file[enrich_file$Cluster == "Tr2",]
barplot(tr2$Score, horiz = T, col = "#FF6666")

tr3 = enrich_file[enrich_file$Cluster == "Tr3",]
barplot(tr3$Score, horiz = T, col = "#FF6666")

tr5 = enrich_file[enrich_file$Cluster == "Tr5",]
barplot(tr5$Score, horiz = T, col = "#FF6666")

# UMAP in tumor clusters (Fig.S4a)
genes_to_show = c("GPC3","ALB","AFP","DLK1","MET","KRT18","KRT19","PROM1","CD44",
                  "VIM","FLT1","CD34","THY1","KIT","OCT4","HMGA2")
FeaturePlot(tumor.integrated, features = genes_to_show, ncol = 5, cols = c("#F7F7F7","#FF0000"))

# dotplot of tumor genes for tumor clusters (Fig.S4b)
DefaultAssay(tumor.integrated) = "RNA"
DotPlot(tumor.integrated, group.by = "integrated_snn_res.0.3", cols = c("#FFFCF5","red","#750906"),
        features = c("GPC3","DLK1","HDAC2","DUSP9","YAP1","CTNNB1","NKD1","COL2A1","EPCAM",
                     "AFP","FOS","JUN","KRT19","SPP1")) +
  theme(axis.title = element_blank(),
        axis.text = element_text(family = "Helvetica", size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1))
