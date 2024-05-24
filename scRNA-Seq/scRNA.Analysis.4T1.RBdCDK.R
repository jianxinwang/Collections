# for Vishnu's single cell data analysis
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

date <- Sys.Date()
project <- '4T1_Dox'

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# human gene to mouse gene mapping
humanMouseMap <- read.csv("/mnt/Jason/RNA-Seq/util/HumanMouseGeneMapping.txt", sep = "\t", header = FALSE)

# convert human gene ids to mouse gene ids for use in cell cycle score calculation
mouse.s.genes <- humanMouseMap[humanMouseMap$V1 %in% cc.genes.updated.2019$s.genes, ]$V2
mouse.g2m.genes <- humanMouseMap[humanMouseMap$V1 %in% cc.genes.updated.2019$g2m.genes, ]$V2

mouse.cd4.genes <- humanMouseMap[humanMouseMap$V1 %in% c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8"), ]$V2

# read input files
count10x.minus.dox <- Read10X_h5("/mnt/Jason/scRNA/data/Yin/COUNT/RS-03901823_4T1-RBdCDKminusDox_RS-03900475_count/outs/filtered_feature_bc_matrix.h5")
count10x.plus.dox <- Read10X_h5("/mnt/Jason/scRNA/data/Yin/COUNT/RS-03901824_4T1-RBdCDKplusDox_RS-03900476_count/outs/filtered_feature_bc_matrix.h5")

sobj.control <- CreateSeuratObject(counts = count10x.minus.dox, project = "Control")
sobj.treated <- CreateSeuratObject(counts = count10x.plus.dox, project = "Treated")

sobjList <- list()
sobjList[['Control']] <- sobj.control
sobjList[['Treated']] <- sobj.treated

sobjFilteredList <- list()

dim(sobj.control)
dim(sobj.treated)

setwd("/mnt/Jason/scRNA/output/Yin")

for (treatment in names(sobjList)){
    sobj <- sobjList[[treatment]]
    
    # calculate % mitocondrial counts
    sobj$percent.mt <- PercentageFeatureSet(sobj, pattern = "^mt-")
    
    # get the counts of the largest expressed genes
    # Malat1 is the largest count gene, so exclude it from counting
    sobj.no.malat <- sobj[rownames(sobj) != 'Malat1', ]
    sobj.no.malat$largest_count <- apply(sobj.no.malat@assays$RNA@counts, 2, max)
    sobj.no.malat$largest_index <- apply(sobj.no.malat@assays$RNA@counts, 2, which.max)
    sobj.no.malat$largest_gene <- rownames(sobj.no.malat)[sobj.no.malat$largest_index]
    sobj.no.malat$percent.largest.gene <- sobj.no.malat$largest_count/sobj.no.malat$nCount_RNA * 100
    
    sobj$largest.gene <- sobj.no.malat$largest_gene
    sobj$percent.largest.gene <- sobj.no.malat$percent.largest.gene
    
    rm(sobj.no.malat)
    
    p <- VlnPlot(sobj, ncol = 4, features=c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.largest.gene")) 
    
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
        width = 5,
        units = 'in',
        res = 300
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
        width = 5,
        units = 'in',
        res = 300
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
        res = 300
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
    
    # get the counts of the largest expressed genes
    # Malat1 is the largest count gene, so exclude it for counting
    sobj.filtered.no.malat <- sobj.filtered[rownames(sobj) != 'Malat1', ]
    sobj.filtered.no.malat$largest_count <- apply(sobj.filtered.no.malat@assays$RNA@counts, 2, max)
    sobj.filtered.no.malat$largest_index <- apply(sobj.filtered.no.malat@assays$RNA@counts, 2, which.max)
    sobj.filtered.no.malat$largest_gene <- rownames(sobj.filtered.no.malat)[sobj.filtered.no.malat$largest_index]
    sobj.filtered.no.malat$percent.largest.gene <- sobj.filtered.no.malat$largest_count/sobj.filtered.no.malat$nCount_RNA * 100
    
    sobj.filtered$largest.gene <- sobj.filtered.no.malat$largest_gene
    sobj.filtered$percent.largest.gene <- sobj.filtered.no.malat$percent.largest.gene
    
    rm(sobj.filtered.no.malat)
    
    # make plots after filtering
    p <- VlnPlot(sobj.filtered, ncol = 4, features=c("nFeature_RNA", "nCount_RNA","percent.mt", "percent.largest.gene")) 
    
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
        width = 5,
        units = 'in',
        res = 300
    )
    print(p)
    dev.off()
    
    sobjFilteredList[[treatment]] <- sobj.filtered
    
    gc()

}   

rm(sobjList)

#####################################################################################################################
# integrated analysis
#####################################################################################################################

# normalize and identify variable features for each dataset independently
sobjFilteredList <- lapply(X = sobjFilteredList, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = sobjFilteredList, nfeatures = 3000)


# identify anchors 
immune.anchors <- FindIntegrationAnchors(object.list = sobjFilteredList, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

# specify that we will perform downstream analysis on the corrected data note that the
# original unmodified data still resides in the 'RNA' assay
DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

# determine how many PCs to include in RunUMAP
p <- ElbowPlot(immune.combined)

png(paste(date, project, "CombinedFilteredElbowPlot.png", sep = '_'),
    width = 5,
    height = 5,
    units = 'in',
    res = 300
)
print(p)
dev.off()

png(paste(date, project, "CombinedFilteredDimHeatmap.png", sep = '_'),
    width = 10,
    height = 7,
    units = 'in',
    res = 300
)
DimHeatmap(immune.combined, dims=1:20, cells=500)
dev.off()

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

names(immune.combined@meta.data)[1] <- "Treatment"


##############################################################################################################
# annotate with ScType
##############################################################################################################

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = immune.combined@assays$integrated@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(immune.combined@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(immune.combined@meta.data[immune.combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(immune.combined@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

immune.combined@meta.data$Annotated = ""
for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    immune.combined@meta.data$Annotated[immune.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# labeled by cluster ID
p3 <- DimPlot(immune.combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

# labeded by annotation
p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'Annotated', pt.size = 0.05) + NoLegend()

p <- p4 + p3


png(paste(date, project, "CombinedFilteredDimPlotUMAP_withAnnotation.png", sep = '_'),
    width = 13,
    height = 7,
    units = 'in',
    res = 300
)
print(p)
dev.off()


# for visualization purpose, we need to use equal number of cells. 

# dim(sobjFilteredList[["Control"]])
# [1] 32285  7338
# > dim(sobjFilteredList[["Treated"]])
# [1] 32285  9419
# > 9419-7338
# [1] 2081

# So we need to down-sample treated to have 7338 cells
set.seed(1234)

sobj.down.sampled <- subset(immune.combined, subset = Treatment == 'Treated')
sobj.down.sampled$barcode <- names(sobj.down.sampled$Treatment)
to.remove <- sample(sobj.down.sampled$barcode, 2081) # 9307-7204
immune.combined$barcode <- names(immune.combined$Treatment)
to.keep <- immune.combined$barcode[!immune.combined$barcode %in% to.remove]

immune.combined.downsampled <- subset(immune.combined, subset = barcode %in% names(to.keep))

rm(sobj.down.sampled)

immune.combined.downsampled.treated <- subset(immune.combined.downsampled, subset = Treatment == 'Treated')
immune.combined.downsampled.control <- subset(immune.combined.downsampled, subset = Treatment == 'Control')
immune.combined.downsampled.treated$Treated <- immune.combined.downsampled.treated$seurat_clusters
immune.combined.downsampled.control$Control <- immune.combined.downsampled.control$seurat_clusters

p1 <- DimPlot(immune.combined.downsampled.control, reduction = "umap", label = TRUE, group.by = "Control") + NoLegend()
p2 <- DimPlot(immune.combined.downsampled.treated, reduction = "umap", label = TRUE, group.by = "Treated") + NoLegend()

p <- p1 + p2

png(paste(date, project, "CombinedFilteredDownSampledDimPlotUMAPByTreatmentSideBySide.png", sep = '_'),
    width = 11,
    height = 7,
    units = 'in',
    res = 300
)
print(p)
dev.off()

# labeled by cluster ID
p3 <- DimPlot(immune.combined.downsampled, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()

# labeded by annotation
p4 <- DimPlot(immune.combined.downsampled, reduction = "umap", label = TRUE, repel = TRUE, 
              group.by = 'Annotated', pt.size = 0.05) + NoLegend()

p <- p4 + p3


png(paste(date, project, "CombinedFilteredDownSampledDimPlotUMAP_withAnnotation.png", sep = '_'),
    width = 11,
    height = 7,
    units = 'in',
    res = 300
)
print(p)
dev.off()


# clean up
rm(immune.combined)

############################################################################################################
# find conserved markers for every cluster compared to all remaining clusters regardless of treatment
############################################################################################################
DefaultAssay(immune.combined.downsampled) <- "RNA"

findConservedMarkers <- function(i) {
    j <- i - 1
    FindConservedMarkers(immune.combined.downsampled, ident.1 = j, grouping.var = "Treatment", verbose = TRUE, max.cells.per.ident = 500)
}

conservedMarkersList <- lapply(1:length(levels(immune.combined.downsampled$seurat_clusters)), findConservedMarkers)

rownames(conservedMarkersList[[1]][1:4, ])

getTopTwoConservedMarkers <- function(i){
  rownames(conservedMarkersList[[i]][1:2, ])
}

topTwoConserverdMarkers <- unlist(lapply(1:length(levels(immune.combined.downsampled$seurat_clusters)), getTopTwoConservedMarkers))
topTwoConserverdMarkers <- unique(topTwoConserverdMarkers)

# show conserved markers in FeaturePlots
p <- FeaturePlot(immune.combined.downsampled, ncol = 8, features = topTwoConserverdMarkers)

png(paste(date, project, "Top2ConservedMarkersFeaturePlot.png", sep = '_'),
    width = 32,
    height = 20,
    units = 'in',
    res = 100
)
print(p)
dev.off()


#
# DotPlot to visulize the results
#
DefaultAssay(immune.combined.downsampled) <- "RNA"
Idents(immune.combined.downsampled) <- "Annotated"
markers.to.plot <- topTwoConserverdMarkers

png(paste(date, project, "ConservedMarersDotPlot.png", sep = '_'),
    width = 14,
    height = 8,
    units = 'in',
    res = 200)

DotPlot(immune.combined.downsampled, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "Treatment") +
  RotatedAxis()
dev.off()

############################################################################################################
# find markers for every cluster compared to all remaining clusters
############################################################################################################

DefaultAssay(immune.combined.downsampled) <- "integrated"
Idents(immune.combined.downsampled) <- immune.combined.downsampled$seurat_clusters
immune.combined.downsampled.markers <- FindAllMarkers(immune.combined.downsampled, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 500)

top10 <- immune.combined.downsampled.markers %>% group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) 

##################################
# Top feature heatmap
##################################

p <- DoHeatmap(immune.combined.downsampled, features = top10$gene) + NoLegend()

png(paste(date, project, "CombinedFilteredTop10FeatureHeatmap.png", sep = '_'),
    width = 17,
    height = 17,
    units = 'in',
    res = 200
)
print(p)
dev.off()

###########################################################################################################################
# Identify differentially expressed markers between treatment conditions. Note: there is no biological replicates. The 
# difference reported could be due to treatment or due to difference in samples even with identical treatments, or both.
###########################################################################################################################


#
# Method 1: scatterplot to show outlier genes
#
DefaultAssay(immune.combined.downsampled) <- "RNA"
immune.combined.downsampled$celltype <- immune.combined.downsampled$Annotated
Idents(immune.combined.downsampled) <- "celltype"

makeScatterplot <- function(x){
  cellType <- gsub(' +', '_', x)
  
  cells.of.interest <- subset(immune.combined.downsampled, idents = x)
  
  Idents(cells.of.interest) <- "Treatment"
  
  avg.exp <- as.data.frame(log1p(AverageExpression(cells.of.interest, verbose = FALSE)$RNA))
  avg.exp$gene <- rownames(avg.exp)
  
  
  # calculate ratio between treatments for outlier labeling
  avg.exp$min <- apply(avg.exp, 1, min)
  avg.exp$ratio <- avg.exp$Treated/avg.exp$Control
  avg.exp <- na.omit(avg.exp)
  avg.exp <- avg.exp[!is.infinite(avg.exp$ratio), ]
  
  down.genes <- avg.exp[avg.exp$min > 0.5 & avg.exp$ratio < 0.8, ]
  down.genes <- down.genes[order(down.genes$ratio, decreasing = TRUE), ]
  
  
  up.genes <- avg.exp[avg.exp$min > 0.5 & avg.exp$ratio > 1.2, ]
  up.genes <- up.genes[order(up.genes$ratio, decreasing = TRUE), ]
  
  genes.to.label = c(up.genes$gene[1:20], down.genes$gene[1:20])
  
  p1 <- ggplot(avg.exp, aes(Control, Treated)) + geom_point(size = 0.5) + ggtitle(x)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0)
  
  png(paste(date, project, cellType, "NormalizedCountScatterPlot.png", sep = '_'),
      width = 6,
      height = 6,
      units = 'in',
      res = 300)
  print(p1)
  dev.off()
  
}

# The last element has less than three cells so skip it
rv <- lapply(unique(immune.combined.downsampled$Annotated), makeScatterplot)
#rv <- lapply(head(unique(immune.combined.downsampled$Annotated), -1), makeScatterplot)

# 
# Method 2: calculate differential genes using FindMarker function
#
immune.combined.downsampled$celltype.condition <- paste(immune.combined.downsampled$Annotated, immune.combined.downsampled$Treatment, sep = "_")
immune.combined.downsampled$celltype <- immune.combined.downsampled$Annotated
Idents(immune.combined.downsampled) <- "celltype.condition"

getDeGeneBetweenConditions2 <- function(x){
  cellType <- gsub(' +', '_', x)
  
  print(cellType)
  
  ident1 <- paste(x, 'Treated', sep = '_')
  ident2 <- paste(x, 'Control', sep = '_')
  
  
  
  response <- FindMarkers(immune.combined.downsampled, ident.1 = ident1, ident.2 = ident2,
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

cellTypes <- unique(immune.combined.downsampled$Annotated)
#cellTypes <- cellTypes[which(cellTypes != "ISG expressing immune cells")]

rv <- lapply(cellTypes, getDeGeneBetweenConditions2)

genes.of.interest <- unlist(rv)
genes.of.interest <- unique(genes.of.interest[!is.na(genes.of.interest)])

####################################################################################
# Violin plots for select de genes
####################################################################################


plots <- VlnPlot(immune.combined.downsampled, features = genes.of.interest, 
                 split.by = "Treatment", group.by = "Annotated",
                 pt.size = 0, combine = FALSE, assay = 'RNA', slot = "data")

p <- wrap_plots(plots = plots, ncol = 7)

png(paste(date, project, "DeGenesViolinPlots.png", sep = '_'),
    width = 35,
    height = 22,
    units = 'in',
    res = 200)
print(p)
dev.off()

####################################################################################
# compare cell count between treatment conditions for each of the clusters
####################################################################################
barcode.control <- names(immune.combined.downsampled$Annotated[immune.combined.downsampled$Treatment == 'Control'])
barcode.treated <- names(immune.combined.downsampled$Annotated[immune.combined.downsampled$Treatment == 'Treated'])

count.control <- table(immune.combined.downsampled$Annotated[names(immune.combined.downsampled$Annotated) %in% barcode.control])
count.treated <- table(immune.combined.downsampled$Annotated[names(immune.combined.downsampled$Annotated) %in% barcode.treated])

count.ratio <- data.frame(Control  = count.control, Treated = count.treated)
count.ratio <- count.ratio[, c(1,2,4)]
colnames(count.ratio) <- c("CellType", "Control", "Treated")

data.for.barplot <- melt(count.ratio, id.vars = "CellType", variable.name = "Treatment", value.name = "Count")

p <- ggplot(data.for.barplot, aes(x = CellType, y = Count, fill = Treatment)) +
    geom_bar(stat="identity", position=position_dodge()) + theme(axis.text.x = element_text(angle=45, vjust=1.0, hjust=1)) +
    theme(legend.position = c(0.7, 0.85))


png(paste(date, project, "CellTypeCountByTreatmentBarplot.png", sep = '_'),
    width = 7,
    height = 5,
    units = 'in',
    res = 200
)
print(p)
dev.off()

count.ratio$`Treated/Control` <- count.ratio$Treated/count.ratio$Control
count.ratio$`Vehile/Treated` <- count.ratio$Control/count.ratio$Treated

#### WRITE SESSION INFO TO FILE ####################################################################
writeLines(capture.output(sessionInfo()), paste(date, "scRNA-Seq_Analysis_SessionInfo.txt", sep = '_'))

