library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)


setwd("/mnt/Jason/scRNA/output/Vishnu")

#
# Step -1: Convert data from Seurat to Python / anndata
#

# load previously saved image
load("/mnt/Jason/scRNA/output/Vishnu/AKB6.Rdata")

immune.combined.downsampled$barcode <- colnames(immune.combined.downsampled)
immune.combined.downsampled$UMAP_1 <- immune.combined.downsampled@reductions$umap@cell.embeddings[,1]
immune.combined.downsampled$UMAP_2 <- immune.combined.downsampled@reductions$umap@cell.embeddings[,2]
write.csv(immune.combined.downsampled@meta.data, file='metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(immune.combined.downsampled, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(immune.combined.downsampled@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write gene names
write.table(
    data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
    quote=F,row.names=F,col.names=F
)
