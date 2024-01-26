library(Seurat)
library(dplyr)
library(ggplot2)

#Loading data (use either cardiac or CKO dataset)
countData <- read.csv("cardiac/output/cardiac_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("cardiac/output/cardiac_sc_celltype.csv", header = TRUE)
cell_type <- cell_type[,2:3]
rownames(cell_type) <- cell_type[,1]
cardiac_cells <- CreateSeuratObject(counts = countData, project = "cardiac_project", min.cells = 3, min.features = 200, meta.data = cell_type)

#Calculate percentage of mitochondrial genes in each cell
cardiac_cells[["percent.mt"]] <- PercentageFeatureSet(cardiac_cells, pattern = "mt-")
VlnPlot(cardiac_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(cardiac_cells, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cardiac_cells, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Remove cells that have a given number of genes OR
#more than a certain % of mitochondrial genes
cardiac_cells <- subset(cardiac_cells, subset = nFeature_RNA > 200 & percent.mt < 25)

#Normalize data with LogNormalize
cardiac_cells <- NormalizeData(cardiac_cells, normalization.method = "LogNormalize", scale.factor = 10000)

#Find highly variable genes 
#(i.e. highly expressed in some cells, and lowly expressed in others)
cardiac_cells <- FindVariableFeatures(cardiac_cells, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cardiac_cells), 10)
top10

# plot variable features with and without labels
plot <- VariableFeaturePlot(cardiac_cells)
plot1 <- LabelPoints(plot = plot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
plot1

#Scale data
all.genes <- rownames(cardiac_cells)
cardiac_cells <- ScaleData(cardiac_cells, features = all.genes)

#Perform linear dimension reduction...
cardiac_cells <- RunPCA(cardiac_cells, features = VariableFeatures(object = cardiac_cells))

#...And visualize PCA in some ways
VizDimLoadings(cardiac_cells, dims = 1:2, reduction = "pca")
DimPlot(cardiac_cells, reduction = "pca")
DimHeatmap(cardiac_cells, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(cardiac_cells, dims = 1:15, cells = 500, balanced = TRUE)

#Determine the 'dimensionality' of the dataset
#'Significant' PCs will show a strong enrichment of 
#features with low p-values (solid curve above the dashed line).
cardiac_cells <- JackStraw(cardiac_cells, num.replicate = 100)
cardiac_cells <- ScoreJackStraw(cardiac_cells, dims = 1:20)

#Visualise with Straw plot
JackStrawPlot(cardiac_cells, dims = 1:20)

#Or with Elbow plot
ElbowPlot(cardiac_cells)

#Cluster cells and look at the first 5 cell types
cardiac_cells <- FindNeighbors(cardiac_cells, dims = 1:7)
cardiac_cells <- FindClusters(cardiac_cells, resolution = 0.1)
head(Idents(cardiac_cells), 8)

#Run UMAP/tSNEs
cardiac_cells <- RunUMAP(cardiac_cells, umap.method = 'umap-learn', metric = 'correlation', dims = 1:7)
cardiac_cells <- RunTSNE(cardiac_cells, dims = 1:7)
umap <- DimPlot(cardiac_cells, reduction = "umap", label = T, label.box = F, repel =T)
tsne <- DimPlot(cardiac_cells, reduction = "tsne", label = T, label.box = F, repel= T)

umap
tsne

# Save Unlabelled clusterings
ggsave(tsne, file = "Unlabelled tSNE V3.png", width = 7.5, height = 6, units = 'in', dpi = 600)

#Save object so far
saveRDS(cardiac_cells, file = "sc_cardiac_obj.rds")

# Find markers for every cluster compared to all remaining cells, report only the positive
# ones
cardiac_cells.markers <- FindAllMarkers(cardiac_cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cardiac_cells.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)%>% 
  write.csv(file='cell_markers_res0.1_dims7.csv')
cluster0.markers <- FindMarkers(cardiac_cells, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

  
FeaturePlot(cardiac_cells, features = c("Acta2", 'Acsl1'), reduction="tsne")

#Get markers of the tree CMs subpopulations
Cms_1 <- FindMarkers(cardiac_cells, ident.1 = 'CMs 1', ident.2 = c("CMs 2", "CMs 3"))
Cms_2 <- FindMarkers(cardiac_cells, ident.1 = 'CMs 2', ident.2 = c("CMs 1", "CMs 3"))
Cms_3 <- FindMarkers(cardiac_cells, ident.1 = 'CMs 3', ident.2 = c("CMs 2", "CMs 1"))

write.csv(Cms_3, file = 'CMs 3.csv')