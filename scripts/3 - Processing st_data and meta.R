library(Seurat)

#Load Spatial reference dataset #1 (NEONATAL MOUSE HEART)
heart_obj <- Read10X("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Spatial REF/Human left ventricle/")
annotations <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Spatial REF/Neonatal Mouse Heart/Mock_Heart_D7PI/spatial/tissue_positions_list.csv")

#Load Spatial reference datatset #2 (10 WK Male mice Sham for MI)
heart_obj <- Read10X("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Spatial REF/10 Week male mice Sham/")
annotations <- read.csv("C:/Users/chelu/OneDrive - The University of Manchester/PhD OneDrive/Coding/Python Material/Bulk2space/Datasets/Spatial REF/10 Week male mice Sham/GSM5355663_WT_Sham_tissue_positions_list.csv")

# Filter for spots covered by tissue and normalize counts
annotations <- subset(annotations, in_tissue == 1)
heart_norm <- NormalizeData(heart_obj,
                            normalization.method = 'LogNormalize')

# Generate st_data
st_data <- heart_norm[,annotations[,1]]
colnames(st_data) <- paste('spot', 1:nrow(annotations['barcode']), sep='_')

# Generate st_meta
st_meta <- as.data.frame(annotations[2:4])
colnames(st_meta)[1] <- 'Spot' 
rownames(st_meta) <- colnames(st_data)
st_meta['Spot'] <- colnames(st_data)

# Write st_data and st_meta
write.csv(st_data, file="transverse_st_data.csv")
write.csv(st_meta, file="transverse_st_meta.csv")

