library(Seurat)
library(dplyr)
library(plyr)
library(ggplot2)
library(plyr)
library(SeuratWrappers)
library(RColorBrewer)
library(tibble)
library(ggpubr)

setwd("/data/aronow/Kang/SingleCellHBPilot-HB17/results/cellxgene/")

### load data
filtered_combined.3types <- readRDS("combined.3types_v2.RDS")
merged.tumors <- readRDS("merged.tumors_v2d.RDS")
integrated_hepatoblastoma <- readRDS("hepatoblastoma_v7.RDS")
brks <- c(0,0.25,0.5,0.75,1)


### Fig.3A Plotting cell class on UMAP
integrated_hepatoblastoma$cell_type <- mapvalues(integrated_hepatoblastoma$cell_type, from=c("NK.T cells"),to=c("NK cells"))
integrated_hepatoblastoma$cell_type <- mapvalues(integrated_hepatoblastoma$cell_type, 
                                                 from=c("Endothelial cells","Hepatocytes","Cholangiocytes","Kuffer cells","Hepatic Stellate cells","T cells","NK cells","B cells","Tumor cells","Inflammatory Monocytes"),
                                                 to=c("2,Endothelial cells","1,Hepatocytes","9,Cholangiocytes","3,Kupffer cells","4,Hepatic Stellate cells","5,T cells","6,NK cells","7,B cells","0,Tumor cells","8,Inflammatory Monocytes"))
colors = c("#FFC0CB","#32CD32","#FFD700","#00BFFF","#FFA500","#BA55D3","#D8BFD8","#CCE5FF","#FA8072","#990000")
DimPlot(integrated_hepatoblastoma,group.by = "cell_type",cols = colors) + 
        theme(legend.text = element_text(size=25,face = "bold",family = "Helvetica"),
              legend.key = element_rect(size=10), legend.key.size = unit(2.0, "lines"))



### Fig.3B Expression of key genes on UMAP
FeaturePlot(merged.tumors, features = c("rna_GPC3")) + theme(plot.title =element_blank()) + 
            scale_color_gradient2(low = "#FFFFFF", mid="#FFFF99", high = "#CC0000", midpoint=9, guide="colourbar")
FeaturePlot(merged.tumors, features = c("rna_ALB")) + theme(plot.title =element_blank()) + 
  scale_color_gradient2(low = "#FFFFFF", mid="#FF9999", high = "#CC0000", midpoint=8, guide="colourbar")


### Fig.3C Violin plot for tumor genes
VlnPlot(filtered_combined.3types,features = c("CTNNB1"), group.by = "cell_type_v2",pt.size = 0,y.max = 15) +
        coord_flip() + NoLegend() + theme(panel.border = element_rect(fill=NA,colour="black"),
                                          plot.title = element_blank(),
                                          axis.text = element_blank(),axis.title = element_blank(),
                                          axis.ticks = element_blank())


### Figure 4E. Scatter plot
superbin <- read.table("./superbinexpr.txt",sep="\t",header=1,row.names = 1) # load superbin file with mean values of each gene for background hepatocytes and tumor cells in Tumor/PDX.
superbin <- rename(superbin, c("background.Hepatocytes.6359"="HB Hepatocyte","human.hepatoplastoma.Tumor.cells.7245"="HBT","pdx.Tumor.cells.8027"="PDX"))
superbin$Name <- rownames(superbin)

top25genes <- read.table("./top25genes_hbt_pdx.txt",sep="\t",header=1)
gene_freq = as.data.frame(table(top25genes$id))
shared_ids =gene_freq[gene_freq$Freq>1,"Var1"]
HBT_ids = c()
PDX_ids = c()
tumor_related_ids = c("GPC3","DLK1","EPCAM","KRT19","YAP1") # show specific tumor genes on scatter plot
background_ids = c("ALB","APOC3","HP")
for (i in 1:dim(top25genes))
{
  gene = top25genes$id[i]
  if (top25genes$Sample[i] == "human hepatoplastoma" & !(gene %in% shared_ids)){HBT_ids=append(HBT_ids,gene)}
  if (top25genes$Sample[i] == "pdx" & !(gene %in% shared_ids)){PDX_ids=append(PDX_ids,gene)}
}
# HBT vs. PDX
ggplot(superbin, aes(x=HBT, y=PDX)) + 
        geom_point(size=1,alpha = 0.4) + 
        geom_point(data=superbin[tumor_related_ids,], aes(x=HBT, y=PDX),color="red",size=2) + 
        geom_text(aes(label=ifelse(Name %in% tumor_related_ids,as.character(Name),''),colour="orange",fontface=2),hjust=-0.2,vjust=0.2,size=6) + 
        xlab(expression("Expression (HBT tumor cells)")) + ylab("Expression (PDX tumor cells)") + 
        theme(panel.background = element_rect(fill = NA, color = "black"),
              axis.title = element_text(size=20, family = "Helvetica"),
              axis.text = element_text(size = 15, family = "Helvetica"),
              legend.position = "none")
lm1 = lm(superbin$PDX ~ superbin$HBT)
summary(lm1)$r.squared
# HBT vs. HBB
ggplot(superbin, aes(x=`HB Hepatocyte`, y=HBT)) + 
  geom_point(size=1,alpha = 0.4) + 
  geom_point(data=superbin[c(tumor_related_ids,background_ids),], aes(x=`HB Hepatocyte`, y=HBT),color="red",size=2) + 
  geom_text(aes(label=ifelse(Name %in% c(tumor_related_ids,background_ids),as.character(Name),''),colour="orange",fontface=2),hjust=0.4,vjust=-0.2,size=6) + 
  xlab(expression("Expression (HBB hepatocytes)")) + ylab("Expression (HBT tumor cells)") + 
  theme(panel.background = element_rect(fill = NA, color = "black"),
        axis.title = element_text(size=20, family = "Helvetica"),
        axis.text = element_text(size = 15, family = "Helvetica"),
        legend.position = "none")
lm2 = lm(superbin$HBT ~ superbin$`HB Hepatocyte`)
summary(lm2)$r.squared
# PDX vs. HBB
ggplot(superbin, aes(x=`HB Hepatocyte`, y=PDX)) + 
  geom_point(size=1,alpha = 0.4) + 
  geom_point(data=superbin[c(tumor_related_ids,background_ids),], aes(x=`HB Hepatocyte`, y=PDX),color="red",size=2) + 
  geom_text(aes(label=ifelse(Name %in% c(tumor_related_ids,background_ids),as.character(Name),''),colour="orange",fontface=2),hjust=0.4,vjust=-0.2,size=6) + 
  xlab(expression("Expression (HBB hepatocytes)")) + ylab("Expression (PDX tumor cells)") + 
  theme(panel.background = element_rect(fill = NA, color = "black"),
        axis.title = element_text(size=20, family = "Helvetica"),
        axis.text = element_text(size = 15, family = "Helvetica"),
        legend.position = "none")
lm3 = lm(superbin$PDX ~ superbin$`HB Hepatocyte`)
summary(lm3)$r.squared



### Fig.5A,B Tumor clusters and sources on umap
DimPlot(merged.tumors, group.by = "cluster",cols = c("#3399FF","#00CC66")) + xlim(-6,8.5)+ 
  theme(axis.text = element_text(size=20,family="Helvetica"),
        axis.title = element_text(size=20,family="Helvetica"),
        legend.position = c(0.7,0.92),legend.text = element_text(size=25,family = "Helvetica"),
        legend.key = element_rect(size=10), legend.key.size = unit(2, "lines"))

DimPlot(merged.tumors, group.by = "tissue",cols = c("#3399FF","#00CC66")) + xlim(-6,8.5)+ 
  theme(axis.text = element_text(size=20,family="Helvetica"),
        axis.title = element_text(size=20,family="Helvetica"),
        legend.position = c(0.7,0.92),legend.text = element_text(size=25,family = "Helvetica"),
        legend.key = element_rect(size=10), legend.key.size = unit(2, "lines"))


### Fig.5B Distribution of sources for each cluster
dist_cells <- as.data.frame(table(merged.tumors$tissue, merged.tumors$seurat_clusters))
colnames(dist_cells) <- c("tissue","cluster","num of cells")
ggplot(dist_cells,
       aes(fill=cluster,y=`num of cells`,x=tissue)) +
  geom_bar(stat = "identity",position="fill") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("tissue") + ylab("fraction of cells") + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, family = "Helvetica",size=20,face="bold"),panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),axis.line = element_line(),axis.title = element_blank(),
        legend.text = element_text(size=25),legend.title = element_text(size=25),
        plot.margin = unit(c(1,1,1,1),"cm")) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks))


### Fig.S1A  umap for background liver and HB T
background <- readRDS("background_v2.RDS")
background <- background[,background$cell_type!="Dead cells"]
HBT <- readRDS("tumor_v2.RDS")

DimPlot(background,group.by = "cell_type") + NoLegend() + 
  theme(axis.title=element_text(size=20),axis.text = element_text(size=20))

DimPlot(HBT,group.by = "cell_type") + NoLegend() + 
  theme(axis.title=element_text(size=20),axis.text = element_text(size=20))


### Fig.S1B Cell distribution of cell types across samples
dist_cells <- as.data.frame(table(integrated_hepatoblastoma$tissue_v2, integrated_hepatoblastoma$cell_type))
colnames(dist_cells) <- c("sample","cell type","num of cells")
dist_cells <- dist_cells[dist_cells$sample != "PDX",]
dist_cells <- dist_cells[dist_cells$`cell type` != "Tumor cells",]
dist_cells$`cell type` <- mapvalues(dist_cells$`cell type`, from=c("NK cells","T cells"),to=c("NK.T cells","NK.T cells"))

ggplot(dist_cells,
       aes(fill=`cell type`,y=`num of cells`,x=sample)) +
  geom_bar(stat = "identity",position="fill") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme(axis.text.x = element_text(angle = 35, hjust = 1, family = "Helvetica",size=20,face="bold"),panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),axis.line = element_line(),axis.title = element_blank(),
        legend.text = element_text(size=25),legend.title = element_text(size=25),
        plot.margin = unit(c(1,1,1,1),"cm")) + 
  scale_y_continuous(expand = c(0,0),limits = c(0,NA),breaks = brks, labels=scales::percent(brks))


### Fig.S3 - distribution of cell counts and gene counts
ctable_counts <- read.table("./ctable_counts.txt",sep="\t",header=1,row.names = 1)
ctable_counts$tumor.cluster <- mapvalues(ctable_counts$tumor.cluster,
                                         from = c("C0","C1","C2","C3","C4","C5","C6","C7","C8","C9"),
                                         to = c("T0","T1","T2","T3","T4","T5","T6","T7","T8","T9"))
ctable_counts <- ctable_counts[ctable_counts["tumor.cluster"]!="",]
ggplot(ctable_counts, aes(x=tumor.cluster, y=n_genes, fill = tumor.cluster)) + 
  geom_violin(width=1.5) + ylim(0,7000) + ylab("n_genes") + 
  stat_summary(fun.data = "mean_sdl", mult=1, geom="pointrange", color="#FFCCCC") + 
  theme(axis.text.x = element_text(family = "Helvetica",size=20),panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=12),axis.line = element_line(),axis.title.x = element_blank(),
        axis.title.y = element_text(family = "Helvetica", size = 20),
        legend.text = element_text(size=25),legend.title = element_text(size=25),
        plot.margin = unit(c(1,1,1,1),"cm")) + NoLegend()

ggplot(ctable_counts, aes(x=tumor.cluster, y=n_counts, fill = tumor.cluster)) + 
      geom_violin(width=1.5) + ylim(0,20000) +
      stat_summary(fun.data = "mean_sdl", mult=1, geom="pointrange", color="#FFCCCC") + 
      theme(axis.text.x = element_text(family = "Helvetica",size=20),panel.background = element_blank(),
        axis.text.y = element_text(family = "Helvetica",size=20),axis.line = element_line(),axis.title = element_blank(),
        legend.text = element_text(size=25),legend.title = element_text(size=25),
        plot.margin = unit(c(1,1,1,1),"cm")) + NoLegend()
