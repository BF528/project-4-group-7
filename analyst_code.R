#install packages
library(R.filesets)
library(Seurat)
library(patchwork)
library(dplyr)

#read seurat file
cells <- loadRDS("1st_programmer_deliver.rds")
cbmc <- RunTSNE(cells, dims = 1:10, method = "FIt-SNE")

#find markers for all clusters
t.markers <- FindAllMarkers(cbmc, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)

#plot markers in heatmap
noor <- t.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(cbmc, features = noor$gene)

#assign celltypes to clusters
cluster <- c("Alpha", "Delta", "Beta", "Alpha", "T Cells", "Ductal", "Gamma", "Beta", "Alpha", "Acinar", "Stellate", "Acinar/Ductal", "Macrophage", "Vascular")
names(cluster) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, cluster)

#make tsne plot of clusters
DimPlot(cbmc, reduction = "tsne", label = TRUE, pt.size = 0.7) + NoLegend()

#write marker genes to csv file
write.csv(t.markers, "gene_markers.csv", row.names = FALSE)
