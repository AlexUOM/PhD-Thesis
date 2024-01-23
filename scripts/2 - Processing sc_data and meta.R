#Load data and annotations (Tabula Muris Heart)
heart.file <- readRDS('Tabula Muris/Heart/facs.normalized.Heart.rds')
annotations <- read.csv('tabula-muris-senis-facs-official-raw-obj__cell-metadata__cleaned_ids.csv')

#Filter for RV, LV and 3 months old mice
RV_LV <- colData(heart.file)$subtissue == 'RV' | colData(heart.file)$subtissue == 'LV'
heart.file <- heart.file[,RV_LV]
age <- colData(heart.file)$age_num == 3
heart.file <- heart.file[,age]

#Exclude valve cells
exclude_valve <- colData(heart.file)$free_annotation != "valve cell"
heart.file <- heart.file[,exclude_valve]

#Generate sc_data
sc_data <- as.data.frame(as.matrix(assay(heart.file)))
rownames(annotations) <- annotations[,1]

annotations <- annotations[colnames(sc_data),]

counter <- 0
cell_vector <- character(0)
for (j in 0:length(colnames(sc_data))){
  cell_vector[counter] <- paste("C", j, sep = "_")
  counter <- counter+1
}

colnames(sc_data) <-  cell_vector

#Generate sc_meta
sc_meta <- as.data.frame(colnames(sc_data))
rownames(sc_meta) <- sc_meta[,1]
colnames(sc_meta)[1] <- "Cell"
cell_type <- annotations[,7]
sc_meta["Cell_type"] <- cell_type

#Write sc_data and sc_meta
write.csv(sc_data, file='AC_sc_data.csv')
write.csv(sc_meta, file= 'AC_sc_meta.csv')
