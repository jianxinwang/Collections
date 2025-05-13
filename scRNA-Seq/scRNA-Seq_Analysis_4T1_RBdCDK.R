# This script processes 4T1 RBdCDK single cell RNA-Seq data. Tumors were treated with vehicle (n = 5)
# or Dox (n = 5) for five days and then digested with Liberase (Sigma, 05401020001) and subjected to
# sequencing at GSR using the 10X Chromium and Illumina NovaSeq platform.
library(Seurat)
library(ggplot2)
library(reshape2)
library(sctransform)
library(tidyverse)
library(patchwork)
library(glmGamPoi)
library(SummarizedExperiment)
library(EnhancedVolcano)
library(scFeatureFilter)
library(openxlsx)
library(HGNChelper)
library(ggpubr)
library(viridis)
library(circlize)
library(randomcoloR)
library(ComplexHeatmap)
library(Azimuth)
library(cowplot)
library(SingleR)
library(celldex)
library(ProjecTILs)

date <- Sys.Date()
project <- '4T1_Dox'

# human gene to mouse gene mapping
humanMouseMap <- read.csv("/mnt/Jason/RNA-Seq/util/HumanMouseGeneMapping.txt", sep = "\t", header = FALSE)

# load gene set preparation function. This is for annotation of seurat clusters using scType
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# read input files
count10x.minus.dox <- Read10X_h5("/mnt/Jason/scRNA-Seq/data/Yin/COUNT/RS-03901823_4T1-RBdCDKminusDox_RS-03900475_count/outs/filtered_feature_bc_matrix.h5")
count10x.plus.dox <- Read10X_h5("/mnt/Jason/scRNA-Seq/data/Yin/COUNT/RS-03901824_4T1-RBdCDKplusDox_RS-03900476_count/outs/filtered_feature_bc_matrix.h5")

sobj.control <- CreateSeuratObject(counts = count10x.minus.dox, project = "Control")
sobj.treated <- CreateSeuratObject(counts = count10x.plus.dox, project = "Treated")

sobjList <- list()
sobjList[['Control']] <- sobj.control
sobjList[['Treated']] <- sobj.treated

sobjFilteredList <- list()

dim(sobj.control)
dim(sobj.treated)

setwd("/mnt/Jason/scRNA-Seq/output/Yin")

#########################################################################################################
# QC
#########################################################################################################
for (treatment in names(sobjList)){
    sobj <- sobjList[[treatment]]
    
    # calculate % mitocondrial counts
    sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = "^mt-")
    
    # get the counts of the most highly expressed genes
    # Malat1 is the largest count gene, so exclude it from counting
    sobj.no.malat <- sobj[rownames(sobj) != 'Malat1', ]
    sobj.no.malat$largest_count <- apply(sobj.no.malat@assays$RNA@layers$counts, 2, max)
    sobj.no.malat$largest_index <- apply(sobj.no.malat@assays$RNA@layers$counts, 2, which.max)
    sobj.no.malat$largest_gene <- rownames(sobj.no.malat)[sobj.no.malat$largest_index]
    sobj.no.malat$percent.largest.gene <- sobj.no.malat$largest_count/sobj.no.malat$nCount_RNA * 100
    
    sobj$largest.gene <- sobj.no.malat$largest_gene
    sobj$percent.largest.gene <- sobj.no.malat$percent.largest.gene
    
    rm(sobj.no.malat)
    
    p <- VlnPlot(sobj, ncol = 4, layer = 'counts', features=c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.largest.gene")) 
    
    # save to a file
    png(paste(date, project, treatment, "QC_VlnPlot_before_filtering.png", sep = "_"),
        height = 5,
        width = 10,
        units = 'in',
        res = 100
    )
    print(p)
    dev.off()
    
    # build a qc data matrix
    qc.metrics <- as_tibble(
        sobj[[]],
        rownames="Cell.Barcode"
    ) 
    
    # plot with log10 scale
    p <- qc.metrics %>%
        arrange(percent.mt) %>%
        ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
        geom_point(size=0.7) + scale_colour_viridis(option = "B") +
        ggtitle(paste("QC metrics for", treatment, "before filtering")) +
        scale_x_log10() + scale_y_log10()
    
    png(paste(date, project, treatment, "feature_QC_plot_before_filtering.png", sep = "_"),
        height = 5,
        width = 6,
        units = 'in',
        res = 100
    )
    print(p)
    dev.off()
    
    # get largest gene list
    largest_gene_list <- qc.metrics %>%
        group_by(largest.gene) %>%
        dplyr::count() %>%
        arrange(desc(n)) 
    
    largest_genes_to_plot <- largest_gene_list %>%
        filter(n>150) %>%
        pull(largest.gene)
    
    p <- qc.metrics %>%
        filter(largest.gene %in% largest_genes_to_plot) %>%
        mutate(largest.gene=factor(largest.gene, levels=largest_genes_to_plot)) %>%
        arrange(largest.gene) %>%
        ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest.gene)) +
        geom_point(size=1)+
        scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(12,"Paired")))
    
    png(paste(date, project, treatment, "largest_gene_before_filtering_scatterplot.png", sep = "_"),
        height = 5,
        width = 6,
        units = 'in',
        res = 100
    )
    print(p)
    dev.off()
    
    p <- qc.metrics %>%
        ggplot(aes(percent.largest.gene)) + 
        geom_histogram(binwidth = 0.7, fill="yellow", colour="black") +
        ggtitle("Distribution of Percentage Largest Gene")
    
    png(paste(date, project, treatment, "distribution_of_largest_gene_before_filtering.png", sep = "_"),
        height = 5,
        width = 5,
        units = 'in',
        res = 100
    )
    print(p)
    dev.off()
    
    #
    # filter the dataset based on threshold obtained above
    #
    minFeature <- 200
    maxFeature <- 4000
    maxPercentMT <- 10
    maxPercentLargestGene <- 20
    
    if (treatment == 'treated'){
        maxPercentMT <- 8
        maxPercentLargestGene <- 15
    } 
    
    sobj.filtered <- subset(sobj, subset = nFeature_RNA > minFeature  & 
                                nFeature_RNA < maxFeature &
                                percent.mt < maxPercentMT & percent.largest.gene < maxPercentLargestGene) 
    
    # calculate % mitocondrial counts
    sobj.filtered$percent.mt <- PercentageFeatureSet(sobj.filtered, pattern = "^mt-")
    
    # get the counts of the most highly expressed genes
    # Malat1 is the largest count gene, so exclude it for counting
    sobj.filtered.no.malat <- sobj.filtered[rownames(sobj) != 'Malat1', ]
    sobj.filtered.no.malat$largest_count <- apply(sobj.filtered.no.malat@assays$RNA@layers$counts, 2, max)
    sobj.filtered.no.malat$largest_index <- apply(sobj.filtered.no.malat@assays$RNA@layers$counts, 2, which.max)
    sobj.filtered.no.malat$largest_gene <- rownames(sobj.filtered.no.malat)[sobj.filtered.no.malat$largest_index]
    sobj.filtered.no.malat$percent.largest.gene <- sobj.filtered.no.malat$largest_count/sobj.filtered.no.malat$nCount_RNA * 100
    
    sobj.filtered$largest.gene <- sobj.filtered.no.malat$largest_gene
    sobj.filtered$percent.largest.gene <- sobj.filtered.no.malat$percent.largest.gene
    
    rm(sobj.filtered.no.malat)
    
    # make plots after filtering
    p <- VlnPlot(sobj.filtered, ncol = 4, layer = 'counts', features=c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.largest.gene")) 
    
    # save to a file
    png(paste(date, project, treatment, "QC_VlnPlot_after_filtering.png", sep = "_"),
        height = 5,
        width = 10,
        units = 'in',
        res = 100
    )
    print(p)
    dev.off()
    
    
    # build a qc data matrix
    qc.metrics <- as_tibble(
        sobj.filtered[[]],
        rownames="Cell.Barcode"
    ) 
    # plot with log10 scale
    p <- qc.metrics %>%
        arrange(percent.mt) %>%
        ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) +
        geom_point(size=0.7) + scale_colour_viridis(option = "B") +
        ggtitle(paste("QC metrics for", treatment, "after filtering")) +
        scale_x_log10() + scale_y_log10()
    
    png(paste(date, project, treatment, "feature_QC_plot_after_filtering.png", sep = "_"),
        height = 5,
        width = 7,
        units = 'in',
        res = 100
    )
    print(p)
    dev.off()
    
    sobjFilteredList[[treatment]] <- sobj.filtered
    
    gc()

}   

# clean up
rm(sobjList)
gc()


# for intergrated analysis, we need to use equal number of cells. 

# dim(sobjFilteredList[["Control"]])
# [1] 32285  7338
# > dim(sobjFilteredList[["Treated"]])
# [1] 32285  9439
# > 9419-7338
# [1] 2081

# So we need to down-sample treated to have 7338 cells
set.seed(1234)

sobj.treated <- sobjFilteredList[['Treated']]
sobj.treated$barcode <- names(sobj.treated$orig.ident)

to.keep <- sample(sobj.treated$barcode, 7338)

sobj.treated.downsampled <- subset(sobj.treated, subset = barcode %in% names(to.keep)) 

dim(sobj.treated.downsampled)
# [1] 32285  7338

# merge two Seurat objects
sobj.integrated <- merge(x = sobj.treated.downsampled, y = sobjFilteredList[['Control']])
sobj.integrated$Treatment <- ifelse(sobj.integrated$orig.ident == 'Treated', 'Dox', 'Vehicle')

# run standard analysis workflow
sobj.integrated <- NormalizeData(sobj.integrated)
sobj.integrated <- FindVariableFeatures(sobj.integrated)
sobj.integrated <- ScaleData(sobj.integrated, features = rownames(sobj.integrated))
sobj.integrated <- RunPCA(sobj.integrated)

# determine how many PCs to include in RunUMAP using albowplot
p <- ElbowPlot(sobj.integrated)

png(paste(date, project, "CombinedFilteredElbowPlot.png", sep = '_'),
    width = 5,
    height = 5,
    units = 'in',
    res = 100
)
print(p)
dev.off()

# also visualize with a heatmap
png(paste(date, project, "CombinedFilteredDimHeatmap.png", sep = '_'),
    width = 10,
    height = 7,
    units = 'in',
    res = 150
)
DimHeatmap(sobj.integrated, dims=1:15, cells=500)
dev.off()

# looks like 15 PCs are a good choice
sobj.integrated <- FindNeighbors(sobj.integrated, dims = 1:15, reduction = "pca")
sobj.integrated <- FindClusters(sobj.integrated, resolution = 2, cluster.name = "unintegrated_clusters")
sobj.integrated <- RunUMAP(sobj.integrated, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")


# rename the orig.ident to 'Treatment'
#names(sobj.integrated@meta.data)[1] <- "Treatment"

# Show clusters before integration
png(paste(date, project, "BeforeIntegrationDimPlot.png", sep = '_'),
    width = 12,
    height = 5,
    units = 'in',
    res = 200
)
DimPlot(sobj.integrated, reduction = "umap.unintegrated", group.by = c("Treatment", 'seurat_clusters'))
dev.off()

sobj.integrated <- IntegrateLayers(object = sobj.integrated, method = CCAIntegration, orig.reduction = "pca",
                                               new.reduction = "integrated.cca",
                                               verbose = FALSE)

# re-join layers after integration
# sobj.integrated[["RNA"]] <- JoinLayers(sobj.integrated[["RNA"]])

sobj.integrated <- FindNeighbors(sobj.integrated, reduction = "integrated.cca", dims = 1:15)
sobj.integrated <- FindClusters(sobj.integrated, resolution = 3.0)
sobj.integrated <- RunUMAP(sobj.integrated, dims = 1:15, reduction = "integrated.cca")

# Show clusters after integration
png(paste(date, project, "AfterIntegrationDimPlot.png", sep = '_'),
    width = 10,
    height = 5,
    units = 'in',
    res = 200
)
DimPlot(sobj.integrated, reduction = "umap", group.by = c("Treatment", 'seurat_clusters'), label = T)
dev.off()


# Correct cell cycle phase caused cell difference
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert to mouse orthologs
s.genes.mouse <- unique(humanMouseMap$V2[humanMouseMap$V1 %in% s.genes])
g2m.genes.mouse <- unique(humanMouseMap$V2[humanMouseMap$V1 %in% g2m.genes])

sobj.integrated <- CellCycleScoring(sobj.integrated, s.features = s.genes.mouse, 
                                    g2m.features = g2m.genes.mouse, set.ident = TRUE)

# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
sobj.integrated <- RunPCA(sobj.integrated, features = c(s.genes.mouse, g2m.genes.mouse))

# Save this DimPlot to show cell cycle phase of the cancer cluster
png(paste(date, project, "CellCyclePhaseDimPlot.png", sep = '_'),
    width = 6,
    height = 5,
    units = 'in',
    res = 200
)
DimPlot(sobj.integrated, reduction = "umap",)
dev.off()


##############################################################################################################
# annotate with ScType
##############################################################################################################

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,
                         # Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

seurat_package_v5 <- isFALSE('counts' %in% names(attributes(sobj.integrated[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(sobj.integrated[["RNA"]]$scale.data) else as.matrix(sobj.integrated[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(sobj.integrated@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(sobj.integrated@meta.data[sobj.integrated@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sobj.integrated@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])

# add annotation
sobj.integrated@meta.data$Annotation = ""

for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  sobj.integrated@meta.data$Annotation[sobj.integrated@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# reorder the levels
sobj.integrated$Treatment <- factor(sobj.integrated$Treatment, levels = c('Vehicle', 'Dox'))
p2 <- DimPlot(sobj.integrated, reduction = "umap", split.by = "Treatment", label = FALSE)

png(paste(date, project, "CombinedFilteredDownSampledDimPlotUMAPByTreatmentSideBySide.png", sep = '_'),
    width = 14,
    height = 7,
    units = 'in',
    res = 150
)
print(p2)
dev.off()

# labeled by cluster ID
p3 <- DimPlot(sobj.integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

# labeled by annotation
p4 <- DimPlot(sobj.integrated, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'Annotation', pt.size = 0.05) #+ NoLegend()

p <- p4 + p3


png(paste(date, project, "CombinedFilteredDownsampledDimPlotUMAP_withAnnotation.png", sep = '_'),
    width = 17.5,
    height = 7,
    units = 'in',
    res = 150
)
print(p)
dev.off()

# try different annotation database
ref.immune <- ImmGenData()

sobj.integrated <- JoinLayers(sobj.integrated)

sobj.integrated.sce <- as.SingleCellExperiment(sobj.integrated)


pred.immune <- SingleR(test = sobj.integrated.sce, ref = ref.immune, assay.type.test=1,
                       labels = ref.immune$label.main)

sobj.integrated@meta.data$Annotation2 = pred.immune$labels[match(rownames(sobj.integrated@meta.data), rownames(pred.immune))]

p5 <- DimPlot(sobj.integrated, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'Annotation2', pt.size = 0.05) # + NoLegend()

p <- p5 + p3
png(paste(date, project, "CombinedFilteredDownsampledDimPlotUMAP_withAnnotation2.png", sep = '_'),
    width = 17,
    height = 7,
    units = 'in',
    res = 150
)
print(p)
dev.off()



# Macrophages cluster have several subclusters, so indentify what they are with FeaturePlot
macrophage.markers.m1 <- c('nos2', 'Il1a', 'Il1b', 'Il6', 'Tlr2', 'Tlr4', 'Cd80', 'Cd86')
macrophage.markers.m2 <- c( 'Trem2', 'Cd115', 'Cd206', 'Pparg', 'Arg1', 'Cd163', 'Cd301', 'Clec7a', 'Pdcd1lg2', 'Retnla')
tam <- c('Ccr2', 'Csf1r', 'Marco', 'Pdcd1lg2', 'Cd40', 'Ccl2', 'Csf1', 'Cd16', 'Pdgfb')


png(paste(date, project, "M1_MacrophageFeaturePlot.png", sep = '_'),
    width =7,
    height = 25,
    units = 'in',
    res = 200
)
FeaturePlot(sobj.integrated, features = macrophage.markers.m1, split.by = 'Treatment',
                 min.cutoff = "q0", max.cutoff = "q60", reduction = 'umap')

dev.off()

png(paste(date, project, "M2_MacrophageFeaturePlot.png", sep = '_'),
    width =7,
    height = 25,
    units = 'in',
    res = 200
)
FeaturePlot(sobj.integrated, features = macrophage.markers.m2, split.by = 'Treatment',
            min.cutoff = "q0", max.cutoff = "q60", reduction = 'umap')

dev.off()

png(paste(date, project, "TAM_MacrophageFeaturePlot.png", sep = '_'),
    width =7,
    height = 25,
    units = 'in',
    res = 200
)
FeaturePlot(sobj.integrated, features = tam, split.by = 'Treatment',
            min.cutoff = "q0", max.cutoff = "q60", reduction = 'umap')

dev.off()


#
# Based on annotations above, some seurat clusters can be renames as below
#

# make a seurat_cluster id to annotation mapping df
seurat.to.annot <- list()
seurat.to.annot[['0']] <- 'Naive CD4+'
seurat.to.annot[['1']] <- 'Neutrophils'
seurat.to.annot[['2']] <- 'Neutrophils'
seurat.to.annot[['3']] <- 'Neutrophils'
seurat.to.annot[['4']] <- 'Neutrophils'
seurat.to.annot[['5']] <- 'ISG expressing'
seurat.to.annot[['6']] <- "Naive B"
seurat.to.annot[['7']] <- "Naive CD8+"
seurat.to.annot[['8']] <- "Myeloid DC"
seurat.to.annot[['9']] <- 'ISG expressing'
seurat.to.annot[['10']] <- 'ISG expressing'
seurat.to.annot[['11']] <- 'Neutrophils'
seurat.to.annot[['12']] <- 'Neutrophils'
seurat.to.annot[['13']] <- 'M2 Macrophages'
seurat.to.annot[['14']] <- 'Memory CD4+'
seurat.to.annot[['15']] <- 'Unknown'
seurat.to.annot[['16']] <- 'Naive B'
seurat.to.annot[['17']] <- 'M2 Macrophages'
seurat.to.annot[['18']] <- 'Cancer'
seurat.to.annot[['19']] <- 'M1 Macrophages'
seurat.to.annot[['20']] <- 'M1 Macrophages'
seurat.to.annot[['21']] <- 'Neutrophils'
seurat.to.annot[['22']] <- 'ISG expressing'
seurat.to.annot[['23']] <- 'M2 Macrophages'
seurat.to.annot[['24']] <- 'TAM'
seurat.to.annot[['25']] <- 'Memory CD4+'
seurat.to.annot[['26']] <- 'Pre-B'
seurat.to.annot[['27']] <- 'Cancer'
seurat.to.annot[['28']] <- 'Tregs'
seurat.to.annot[['29']] <- 'M2 Macrophages'
seurat.to.annot[['30']] <- 'Basophils'
seurat.to.annot[['31']] <- 'Neutrophils'
seurat.to.annot[['32']] <- 'NK'
seurat.to.annot[['33']] <- 'Endothelial'
seurat.to.annot[['34']] <- 'Myeloid DC'
seurat.to.annot[['35']] <- 'CD4+ KNT-like'
seurat.to.annot[['36']] <- 'Myeloid DC'
seurat.to.annot[['37']] <- 'Memory CD8+'
seurat.to.annot[['38']] <- 'Naive B'
seurat.to.annot[['39']] <- 'Pre-B'
seurat.to.annot[['40']] <- "Naive CD4+"

Idents(sobj.integrated) <- 'seurat_clusters'
new.cluster.ids <- paste(unlist(seurat.to.annot))
names(new.cluster.ids) <- names(seurat.to.annot)
sobj.integrated <- RenameIdents(sobj.integrated, new.cluster.ids)

Idents(sobj.integrated) <- factor(Idents(sobj.integrated), 
                                  levels = sort(levels(Idents(sobj.integrated)))[c(1:8, 10:12, 9, 13:19)])
sobj.integrated$Annotation3 <- Idents(sobj.integrated)

p6 <- DimPlot(sobj.integrated, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'Annotation3', pt.size = 0.05)

p <- p6 + p3
png(paste(date, project, "CombinedFilteredDownsampledDimPlotUMAP_ManualAnnotation.png", sep = '_'),
    width = 17,
    height = 7,
    units = 'in',
    res = 150
)
print(p)
dev.off()

# Featureplot for genes specified by Yin
genes.of.interest1 <- c('Cd8a', 'Cd4')
genes.of.interest2 <- c('Foxp3', 'Il2ra')
genes.of.interest3 <- c('Cd44', 'Cd69', 'Icos', 'Prf1', 'Ifng', 'Pdcd1', 'Ctla4')

#Idents(sobj.filtered) <- ''
p <- FeaturePlot(sobj.integrated, features = genes.of.interest1, split.by = 'Treatment',
                 min.cutoff = "q0", max.cutoff = "q80", reduction = 'umap')

png(paste(date, project, "SelectMarkersFeaturePlot2.png", sep = '_'),
    width = 8,
    height = 7.5,
    units = 'in',
    res = 200
)
print(p)
dev.off()

p <- FeaturePlot(sobj.integrated, features = genes.of.interest2, split.by = 'Treatment',
                 min.cutoff = "q0", max.cutoff = "q60", reduction = 'umap')

png(paste(date, project, "SelectMarkersFeaturePlot3.png", sep = '_'),
    width = 8,
    height = 7.5,
    units = 'in',
    res = 200
)
print(p)
dev.off()

p <- FeaturePlot(sobj.integrated, features = genes.of.interest3, split.by = 'Treatment', 
                 min.cutoff = "q0", max.cutoff = "q60", reduction = 'umap')

png(paste(date, project, "SelectMarkersFeaturePlot4.png", sep = '_'),
    width = 7.5,
    height = 25,
    units = 'in',
    res = 100
)
print(p)
dev.off()


#
# DotPlot to visulize the results
#
markers.to.plot <- c(genes.of.interest1, genes.of.interest2, genes.of.interest3)

png(paste(date, project, "SelectMarkersDotPlot.png", sep = '_'),
    width = 8,
    height = 10,
    units = 'in',
    res = 200)

DotPlot(sobj.integrated, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, scale = FALSE, split.by = "Treatment") + RotatedAxis()
dev.off()


############################################################################################################
# find markers for every cluster compared to all remaining clusters
############################################################################################################

sobj.integrated.markers <- FindAllMarkers(sobj.integrated, only.pos = TRUE, 
                                                      min.pct = 0.25, logfc.threshold = 0.25, 
                                                      max.cells.per.ident = 500)

# save to file
write.table(x = sobj.integrated.markers[, 1:6],
            file = paste(date, project, "DifferentialGenesForEachCellType.csv", sep = '_'),
            sep = ',',
            row.names = TRUE,
            col.names = NA
            )

top10 <- sobj.integrated.markers %>% group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) 

######################################################################################################
# Top feature heatmap
######################################################################################################
#Idents(sobj.integrated) <- 'seurat_clusters'
p <- DoHeatmap(sobj.integrated, features = top10$gene)

png(paste(date, project, "CombinedFilteredTop10FeatureHeatmap.png", sep = '_'),
    width = 34,
    height = 30,
    units = 'in',
    res = 100
)
print(p)
dev.off()

###########################################################################################################################
# Identify differentialy expressed markers between treatment conditions
###########################################################################################################################

# 
# calculate differential genes using FindMarker function
#
sobj.integrated$celltype.condition <- paste(sobj.integrated$Annotation3, sobj.integrated$Treatment, sep = "_")
sobj.integrated$celltype <- sobj.integrated$Annotation3
Idents(sobj.integrated) <- "celltype.condition"

getDeGeneBetweenConditions <- function(x){
  cellType <- gsub(' +', '_', x)
  
  print(cellType)
  
  ident1 <- paste(x, 'Dox', sep = '_')
  ident2 <- paste(x, 'Vehicle', sep = '_')
  
  
  
  response <- FindMarkers(sobj.integrated, ident.1 = ident1, ident.2 = ident2,
                          verbose = FALSE, assay = 'RNA')
  
  response <- response[!is.infinite(response$avg_log2FC), ]
  
  # volcano plot
  suppressWarnings(p <- EnhancedVolcano(response,
                                        lab = rownames(response),
                                        x = 'avg_log2FC',
                                        y = 'p_val_adj',
                                        title = "",
                                        subtitle = '',
                                        caption = '',
                                        captionLabSize = 14,
                                        pCutoff = 5e-02,
                                        FCcutoff = 0.25,
                                        pointSize = 1.0,
                                        labSize = 2.5,
                                        colAlpha = 0.5,
                                        legendPosition = 'top',
                                        legendLabSize = 10,
                                        legendIconSize = 4.0,
                                        drawConnectors = TRUE,
                                        widthConnectors = 0.25,
                                        max.overlaps = 20,
                                        arrowheads = FALSE
  ))
  
  png(paste(date, project, cellType, "Treated_vs_Control_Volcanoplot.png", sep = '_'),
      width = 8,
      height = 7,
      units = 'in',
      res = 300
  )
  print(p)
  dev.off()
  
  # convert to Human gene symbols
  tmp <- merge(humanMouseMap, response, by.x = "V2", by.y = 0)
  write.table(x = tmp,
              file = paste(date, project, cellType, "DEGenesBetweenConditions.csv", sep = '_'),
              sep = ',',
              quote = FALSE,
              col.names = NA
  )
  sigGenesUp <- response %>% filter(avg_log2FC > 0.5, p_val_adj < 0.1)
  sigGenesDown <- response %>% filter(avg_log2FC < -0.5, p_val_adj < 0.1)
  
  write.table(x = sigGenesUp,
              file = paste(date, project, cellType, "UpGenesBetweenConditionsFC05Q01.csv", sep = '_'),
              sep = ',',
              quote = FALSE,
              col.names = NA
  )
  write.table(x = sigGenesDown,
              file = paste(date, project, cellType, "DownGenesBetweenConditionsFC05Q01.csv", sep = '_'),
              sep = ',',
              quote = FALSE,
              col.names = NA
  )
  
  genes.of.interest <- c(rownames(sigGenesUp)[1:2], rownames(sigGenesDown)[1:2])
  return(genes.of.interest)
}

cellTypes <- unique(sobj.integrated$Annotation3)

rv <- lapply(cellTypes, getDeGeneBetweenConditions)

####################################################################################
# Violin plots for select de genes
####################################################################################

#sobj.integrated.subset <- subset(sobj.integrated, subset = seurat_clusters %in% c(2, 6, 9, 11, 16, 18))

plots <- VlnPlot(sobj.integrated, features = markers.to.plot, stack = T, flip = T,
                 split.by = "Treatment", group.by = "Annotation3", cols = c("#068993", "#F9776D"),
                 pt.size = 0, combine = FALSE, assay = 'RNA', slot = "data")

p <- wrap_plots(plots = plots, ncol = 1)

png(paste(date, project, "SelectGenesViolinPlots.png", sep = '_'),
    width = 10,
    height = 7,
    units = 'in',
    res = 200)
print(p)
dev.off()

#######################################################################################################
# compare cell count between treatment conditions for each of the clusters
#######################################################################################################
barcode.control <- names(sobj.integrated$Annotation3[sobj.integrated$Treatment == 'Vehicle'])
barcode.treated <- names(sobj.integrated$Annotation3[sobj.integrated$Treatment == 'Dox'])

count.control <- table(sobj.integrated$Annotation3[names(sobj.integrated$Annotation3) %in% barcode.control])
count.treated <- table(sobj.integrated$Annotation3[names(sobj.integrated$Annotation3) %in% barcode.treated])

countData <- data.frame(Control  = count.control, Treated = count.treated)
countData <- countData[, c(1,2,4)]
colnames(countData) <- c("Cell Type", "Vehicle", "Dox")

data.for.barplot <- melt(countData, id.vars = "Cell Type", 
                         variable.name = "Treatment", value.name = "Cell Count")
group.colors <- c(Vehicle = "#333BFF", Dox = "#CC6600")

p <- ggplot(data.for.barplot, aes(x = `Cell Type`, y = `Cell Count`, fill = Treatment)) +
    geom_bar(stat="identity", position=position_dodge()) + theme_classic() + 
    theme(axis.text.x = element_text(angle=45, vjust=1.0, hjust=1)) +
    theme(legend.position.inside = c(0.88, 0.85)) + scale_fill_manual(values=group.colors)
    


png(paste(date, project, "CellTypeCountByTreatmentBarplot.png", sep = '_'),
    width = 5,
    height = 3.5,
    units = 'in',
    res = 200
)
print(p)
dev.off()

# Save data to file
write.table(x = data.for.barplot,
            file = paste(date, project, 'CellTypeCountByTreatment.txt', sep = '_'),
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = T
            )

save(sobj.integrated, file = paste(date, project, "4T1.sobj.integrated.Rdata", sep = '_'))

load("2025-05-07_4T1_Dox_4T1.sobj.integrated.Rdata")

#### WRITE SESSION INFO TO FILE ####################################################################
writeLines(capture.output(sessionInfo()), paste(date, "scRNA-Seq_Analysis_SessionInfo.txt", sep = '_'))



