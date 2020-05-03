install.packages('BiocManager')
BiocManager::install('multtest')
BiocManager::install("hgu95av2.db")
library(biomaRt)
library(Seurat)
library(stringr)
library(dplyr)

###test on example umi matrix#####
#txi <- readRDS("expression_matrix")
##################################
count_matrix <- read.csv("/projectnb/bf528/users/group7/project4/Data/fastq_counts/alevin_matrix/UMI_counts_matrix.csv",header = T)
matrix_rown <- count_matrix[,1]
count_matrix <- count_matrix[,-1]
rownames(count_matrix) <- matrix_rown

##convert ENSG to GENE###
gene_ids_version <- rownames(count_matrix)
gene_ids <- str_replace(gene_ids_version,
                        pattern = ".[0-9]+$",
                        replacement = "")
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl", 
                  host="uswest.ensembl.org",
                  ensemblRedirect = FALSE)
annotLookup <- getBM(
  mart=ensembl,
  attributes=c("ensembl_gene_id", "external_gene_name"),
  filter="ensembl_gene_id",
  values=gene_ids)

matched_gene <- match(gene_ids, annotLookup$ensembl_gene_id)
matched_gene <- annotLookup$external_gene_name[matched_gene]

count_matrix$gene_name <- matched_gene
non_redundent_matrix <- count_matrix %>% group_by(gene_name) %>% summarise_all(funs(sum))
non_redundent_matrix <-na.omit(non_redundent_matrix)
rname <- non_redundent_matrix$gene_name
rownames(non_redundent_matrix) <- rname
non_redundent_matrix$gene_name <- NULL
######
#rownames(count_matrix) <- make.names(matched_gene,unique = T)


pbmc <- CreateSeuratObject(counts = non_redundent_matrix , min.cells = 3, min.features = 200, project = "sc_proj4")

##### violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
LabelPoints(plot = plot1, points = top10, repel = TRUE)



##scale data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


###PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#determine cutoff for dims in PCA
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.6)

#####UMAP
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

pbmc <- RunTSNE(pbmc, dims = 1:10, method = "FIt-SNE")
DimPlot(pbmc, reduction = "tsne")

ct_each_cluster <- table(Idents(pbmc))
ct_each_cluster <- as.data.frame(ct_each_cluster)
ct_each_cluster$percent <- round(100*ct_each_cluster$Freq/sum(ct_each_cluster$Freq), digits = 1)
ct_each_cluster$label = paste(ct_each_cluster$Var1," (", ct_each_cluster$percent,"%)", sep = "")
ct_each_cluster$label2 = paste(ct_each_cluster$Var1," (", ct_each_cluster$Freq,")", sep = "")
cols <- rainbow(nrow(ct_each_cluster))
pie(ct_each_cluster$Freq,labels = ct_each_cluster$label2,col=cols)


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
test_cluster <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
#heatmap generation
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
saveRDS(pbmc,file = "1st_upto_umap_pbmc.rds")


output_test <- readRDS("/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")
