library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

# setwd("/mnt/Jason/scRNA/data/Vishnu")
setwd("/mnt/Jason/scRNA/data/Jin")

# load previously saved image
# load("AKB6.immune.combined.downsampled.Rdata")
load("immune.combined.downsampled.Rdata")


#
# Save data for RNA velocity analysis. Do this separately for vehicle and control samples
#
seurat.obj.control <- subset(immune.combined.downsampled, subset = Treatment == 'Vehicle')
seurat.obj.treated <- subset(immune.combined.downsampled, subset = Treatment == 'Treated')

#
seurat.obj.control$barcode <- colnames(seurat.obj.control)
seurat.obj.control$UMAP_1 <- seurat.obj.control@reductions$umap@cell.embeddings[,1]
seurat.obj.control$UMAP_2 <- seurat.obj.control@reductions$umap@cell.embeddings[,2]
write.csv(seurat.obj.control@meta.data, file='metadata.control.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat.obj.control, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.control.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat.obj.control@reductions$pca@cell.embeddings, file='pca.control.csv', quote=F, row.names=F)

# write gene names
write.table(
    data.frame('gene'=rownames(counts_matrix)),file='gene_names.control.csv',
    quote=F,row.names=F,col.names=F
)


seurat.obj.treated$barcode <- colnames(seurat.obj.treated)
seurat.obj.treated$UMAP_1 <- seurat.obj.treated@reductions$umap@cell.embeddings[,1]
seurat.obj.treated$UMAP_2 <- seurat.obj.treated@reductions$umap@cell.embeddings[,2]
write.csv(seurat.obj.treated@meta.data, file='metadata.treated.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat.obj.treated, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.treated.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat.obj.treated@reductions$pca@cell.embeddings, file='pca.treated.csv', quote=F, row.names=F)

# write gene names
write.table(
    data.frame('gene'=rownames(counts_matrix)),file='gene_names.treated.csv',
    quote=F,row.names=F,col.names=F
)


# Combined
immune.combined.downsampled$barcode <- colnames(immune.combined.downsampled)
immune.combined.downsampled$UMAP_1 <- immune.combined.downsampled@reductions$umap@cell.embeddings[,1]
immune.combined.downsampled$UMAP_2 <- immune.combined.downsampled@reductions$umap@cell.embeddings[,2]
write.csv(immune.combined.downsampled@meta.data, file='metadata.combined.downsampled.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(immune.combined.downsampled, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.combined.downsampled.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(immune.combined.downsampled@reductions$pca@cell.embeddings, file='pca.combined.downsampled.csv', quote=F, row.names=F)

# write gene names
write.table(
    data.frame('gene'=rownames(counts_matrix)),file='gene_names.combined.downsampled.csv',
    quote=F,row.names=F,col.names=F
)

