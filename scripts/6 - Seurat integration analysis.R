library(Seurat)
library(patchwork)
library(metap)
library(ggplot2)

# Load data
# Control
countData <- read.csv("Flox/output/Flox_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("Flox/output/Flox_sc_celltype.csv", header = TRUE)
cell_type <- cell_type[,2:3]
rownames(cell_type) <- cell_type[,1]
control_cells <- CreateSeuratObject(counts = countData, project = "SC_project", min.cells = 3, min.features = 200, meta.data = cell_type)
control_cells <- AddMetaData(control_cells, metadata = "Control", col.name = "Treatment")

# Treatment
countData <- read.csv("CKO/output/CKO_sc_data.csv", header = TRUE, row.names =  1)
cell_type <- read.csv("CKO/output/CKO_sc_celltype.csv", header = TRUE)
cell_type <- cell_type[,2:3]
rownames(cell_type) <- cell_type[,1]
treatment_cells <- CreateSeuratObject(counts = countData, project = "SC_project", min.cells = 3, min.features = 200, meta.data = cell_type)
treatment_cells <- AddMetaData(treatment_cells, metadata = "CKO", col.name = "Treatment")

# Merging the two datasets
datasets.list <- list(control_cells, treatment_cells)

# normalize and identify variable features for each dataset independently
datasets.list <- lapply(X = datasets.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = datasets.list)

# We then identify anchors using the FindIntegrationAnchors() function
# which takes a list of Seurat objects as input
# and use these anchors to integrate the two datasets together with IntegrateData().

cardiac.anchors <- FindIntegrationAnchors(object.list = datasets.list, anchor.features = features)

# this command creates an 'integrated' data assay
cardiac.combined <- IntegrateData(anchorset = cardiac.anchors)

# specify that we will perform downstream analysis on the corrected data, note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(cardiac.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
cardiac.combined <- ScaleData(cardiac.combined, verbose = FALSE)
cardiac.combined <- RunPCA(cardiac.combined)

#Determine the 'dimensionality' of the dataset
#'Significant' PCs will show a strong enrichment of 
#features with low p-values (solid curve above the dashed line).
cardiac.combined <- JackStraw(cardiac.combined, num.replicate = 100)
cardiac.combined <- ScoreJackStraw(cardiac.combined, dims = 1:20)

#Visualise with Straw plot and Elbow plot
JackStrawPlot(cardiac.combined, dims = 1:20)
ElbowPlot(cardiac.combined, ndims = 20)

#Clustering 
cardiac.combined <- FindNeighbors(cardiac.combined, reduction = "pca", dims = 1:7)
cardiac.combined <- FindClusters(cardiac.combined, resolution = 0.01)

# Visualization using UMAP
cardiac.combined <- RunUMAP(cardiac.combined, reduction = "pca", dims = 1:7)
p1 <- DimPlot(cardiac.combined, reduction = "umap", group.by = "Treatment", label.size = 5)
p2 <- DimPlot(cardiac.combined, reduction = "umap", label = TRUE, repel = TRUE, label.size = 5)
p1 + p2
umap <- p1+p2
ggsave(umap, file = "UMAP 2 plots.png", width = 13, height = 5, units = 'in', dpi = 600 )

# Visualization using TSNE
cardiac.combined <- RunTSNE(cardiac.combined, reduction = "pca", dims = 1:7)
p3 <- DimPlot(cardiac.combined, reduction = "tsne", group.by = "Treatment", label.size = 5)
p4 <- DimPlot(cardiac.combined, reduction = "tsne", label = TRUE, repel = TRUE, label.size = 5)
tsne <- p3 + p4
ggsave(tsne, file = "tSNE plots unlaballed.png", width = 13, height = 5, units = 'in', dpi = 600 )

# Remove mini-cluster 7-9 as they might be artifacts
cardiac.combined <- subset(cardiac.combined, idents = c(0,1,2,3,4,5,6))

# Save the seurat object so far
saveRDS(cardiac.combined, file = 'cardiac.combined.rds')

# To visualize the two conditions side-by-side, 
# we can use the split.by argument to show each condition colored by cluster.
side.by.side <- DimPlot(cardiac.combined , reduction = "tsne", split.by = 'Treatment', label.size = 5)
side.by.side
ggsave(side.by.side, file = "tSNE CKO vs Control.png", width = 9, height = 5, units = 'in', dpi = 600, bg= "white" )

#---------------------------------------------------------------------------------------------------------------
# For performing differential expression after integration, we switch back to the original
# data
DefaultAssay(cardiac.combined) <- "RNA"

# DGE is done with FindConservedMarkers() for a single cluster at a time
# Repeat this in a loop for every cluster and generate as many files as the number
# of clusters

cell_types <- levels(cardiac.combined)
for (cell in cell_types){
  markers <- FindConservedMarkers(cardiac.combined, ident.1 = cell, grouping.var = "Treatment")
  cluster_name <- paste(cell, "Cluster", sep= "_")
  file_name <- paste(cluster_name, "csv", sep= ".")
  write.csv(markers, file_name)
}

# We can explore these marker genes for each cluster and use them 
# to annotate our clusters as specific cell types

FeaturePlot(cardiac.combined, reduction = 'tsne',
            features = c("Kcnj8"), split.by = "Treatment")

# Annotate/rename clusters with cell names
cardiac.combined <- RenameIdents(cardiac.combined, `0` = "CMs 1", `1` = "Fibroblasts", 
                                 `2` = "CMs 2", `3` = "Leukocytes", `4` = "SMCs", 
                                `5` = "Endothelial cells", `6` = "Myofibroblasts")
cardiac.combined <- RenameIdents(cardiac.combined, `CMs 1` = '0', `Fibroblasts` = '1', 
                                 `CMs 2` = '2', `Leukocytes` = '3', `SMCs` = '4', 
                                 `Endothelial cells` = '5', `Myofibroblasts` = '6')

# Extract cell proportions
cells <- table(cardiac.combined@meta.data$seurat_clusters, cardiac.combined@meta.data$Treatment)
cells
write.csv(cells, 'Cell proportions.csv')

# Visualise the updated tSNE or UMAP with cell annotations
DimPlot(cardiac.combined, reduction = "tsne", label = TRUE, repel = TRUE)
cardiac.combined <- subset(cardiac.combined, idents = c('CMs 1',
                                                      'Fibroblasts',
                                                      'CMs 2',
                                                      'Leukocytes',
                                                      'SMCs',
                                                      ))

# The DotPlot() function with the split.by parameter can be useful for 
# viewing conserved cell type markers across conditions, showing both the 
# expression level and the percentage of cells in a cluster expressing any given gene. 
# Here we plot 2-3 strong marker genes for each of our 14 clusters.
Idents(cardiac.combined)
Idents(cardiac.combined) <- factor(Idents(cardiac.combined), levels = c("CMs 1", 'CMs 2',"Fibroblasts", 
                                                                        "SMCs", 'Leukocytes', 
                                                                        "Endothelial cells",
                                                                        "Myofibroblasts"      ))
markers.to.plot <- c("Actc1", "Myh6", 'Tnni3', 
                     "Col1a1", "Col3a1", "Gsn",
                     "Acta2", "Myl9", "Tagln",
                     "Lyz2", "C1qa", "Cd68", 
                     "Pecam1", "Vwf", 'Cdh5',
                     "Steap4", "Sdc1", 'Kcnj8')

dotplot <- DotPlot(cardiac.combined, features = markers.to.plot, cols=c('blue', 'red'), dot.scale = 8, col.min = 0, split.by = 'Treatment') +
  RotatedAxis()+
  theme(axis.text.x = element_text(size = 14), # Set the font size of the cell labels
        axis.text.y = element_text(size = 16)) # Set the font size of the gene labels
dotplot
# Saving the plot with good resolution
ggsave(dotplot, file = "Markers dotplot.png", width = 10, height = 8, units = 'in', dpi = 600, bg= 'white' )

#--------------------------------------------------------------------------------------
# Identify differential expressed genes across conditions
library(cowplot)
theme_set(theme_cowplot())

# Average expression of both CKO and Control CMs 1 (Example)
cm.cells <- subset(cardiac.combined, idents = "CMs 1")
Idents(cm.cells) <- "Treatment"
avg.cm.cells <- as.data.frame(log1p(AverageExpression(cm.cells, verbose = FALSE)$RNA))
avg.cm.cells$gene <- rownames(avg.cm.cells)

# Average expression of both CKO and Control fibroblasts (Example)
fibro.cells <- subset(cardiac.combined, idents = "Fibroblasts")
Idents(fibro.cells) <- "Treatment"
avg.fibro.cells <- as.data.frame(log1p(AverageExpression(fibro.cells, verbose = FALSE)$RNA))
avg.fibro.cells$gene <- rownames(avg.fibro.cells)

genes.to.label = c("Ppp1cc",
"Bcap31", "Actg1",
"Vps29",
"Uchl3",
'Jak1',
'Stx12',
'Aup1')
p1 <- ggplot(avg.cm.cells, aes(Control, CKO)) + geom_point() + ggtitle("Cardiomyocytes")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.fibro.cells, aes(Control, CKO)) + geom_point() + ggtitle("Fibroblasts")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2


# We can ask what genes change in different conditions for cells of the same type
# For example, we can find the genes that are different between CKO and control Cardiomyocytes.
# This can be repeated for ALL cell types!
cardiac.combined$celltype.Treatment <- paste(Idents(cardiac.combined), cardiac.combined$Treatment, sep = "_")
cardiac.combined$celltype <- Idents(cardiac.combined)
Idents(cardiac.combined) <- "celltype.Treatment"
for (cell in cell_types){
  CKO_cell <- paste(cell, 'CKO', sep= "_")
  control_cell <- paste(cell, 'Control', sep = "_")
  knockout.response <- FindMarkers(cardiac.combined, 
                                   ident.1 = CKO_cell, 
                                   ident.2 = control_cell)
  knockout.response <- knockout.response[order(knockout.response$avg_log2FC, decreasing = TRUE),]
  file_name <- paste(cell, 'comparison.csv', sep = "_")
  write.csv(knockout.response, file_name)
}

head(knockout.response, n = 15)

# Gene expression can be visualized with FeaturePlot()...
plots <- FeaturePlot(cardiac.combined, reduction = 'tsne', features = c("Xirp2", "Snta1"),
            cols = c('grey', 'red'),
            split.by = 'Treatment')
plots
ggsave(plots, file = "Xirp2 and Snta1 FeaturePlot.png", width = 13, height = 10, units = 'in', dpi = 600, bg= "white" )

#...or VlnPlot()
plots <- VlnPlot(cardiac.combined, features = c('Xirp2', 'Snta1'), split.by = "Treatment", group.by = "celltype",
                 pt.size = 0, combine = FALSE)
violins <- wrap_plots(plots = plots, ncol = 1)
violins
ggsave(violins, file = "Xirp2 and Snta1 Violin.png", width = 13, height = 7, units = 'in', dpi = 600, bg= "white" )

