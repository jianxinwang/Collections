library(pheatmap)
library(dplyr)
library(Seurat)
library(patchwork)
library(SCORPIUS)
library(ggplot2)
library(reshape2)
library(sctransform)
library(tidyverse)


smoothByWindow <- function(x){
    res <- NULL
    y <- c(x[(length(x) - 50):length(x)], x, x[1:50])
    
    for (i in 1:(length(x))){
        new.i <- mean(y[i:(i + 100)])
        res <- c(res, new.i)
    }
    res
}

expDataMCF7_T47D <- Read10X_h5("/mnt/scRNA/data/RS-03616030_ERpos-GEX_RS-03616038_ERpos-FB/filtered_feature_bc_matrix.h5")
expDataHCC1806_MiaPaCa <- Read10X_h5(("/mnt/scRNA/data/RS-03616031_MiaHCC-GEX_RS-03616039_MiaHCC-FB/filtered_feature_bc_matrix.h5"))

dimnames(expDataMCF7_T47D$`Antibody Capture`)[[1]] <- c("MCF7RBposDMSO",
                                                        "MCF7RBposPalbo",
                                                        "MCF7RBneg",
                                                        "T47DRBposDMSO",
                                                        "T47DRBposPalbo",
                                                        "T47DRBneg")

dimnames(expDataHCC1806_MiaPaCa$`Antibody Capture`)[[1]] <- c("HCC1806DMSO",
                                                              "HCC1806Palbo",
                                                              "MiaPaCasgctrlDMSO",
                                                              "MiaPaCasgctrlPalbMRTX",
                                                              "MiaPaCasgRBDMSO",
                                                              "MiaPaCasgRBPalbMRTX")


sobj1 <- CreateSeuratObject(counts = expDataMCF7_T47D$`Gene Expression`, project = "CellCycle") 
sobj1[["HTO"]] <- CreateAssayObject(counts = expDataMCF7_T47D$`Antibody Capture`) 
sobj1 <- NormalizeData(sobj1, assay = "HTO", normalization.method = "CLR")
sobj1 <- HTODemux(sobj1, assay = "HTO", positive.quantile = 0.99)

sobj2 <- CreateSeuratObject(counts = expDataHCC1806_MiaPaCa$`Gene Expression`, project = "CellCycle") 
sobj2[["HTO"]] <- CreateAssayObject(counts = expDataHCC1806_MiaPaCa$`Antibody Capture`) 
sobj2 <- NormalizeData(sobj2, assay = "HTO", normalization.method = "CLR")

sobj1 <- HTODemux(sobj1, assay = "HTO", positive.quantile = 0.99)
sobj2 <- HTODemux(sobj2, assay = "HTO", positive.quantile = 0.99)

# Clean up
# First, we will remove negative cells from the object
sobj1 <- subset(sobj1, idents = "Negative", invert = TRUE)
sobj2 <- subset(sobj2, idents = "Negative", invert = TRUE)

# remove doublet cells
sobj1 <- subset(sobj1, idents = "Doublet", invert = TRUE)
sobj2 <- subset(sobj2, idents = "Doublet", invert = TRUE)


sobjList <- list()
sobjList[['ER']] <- sobj1
sobjList[['MiaHCC']] <- sobj2

dim(sobj1)
#[1] 36601 12033
dim(sobj2)

# loop through each hashTag and performs downstream analysis
setwd("/mnt/scRNA/output")

for (celllines in names(sobjList)){
    sobj <- sobjList[[celllines]]
    
    for (treatment in unique(sobj$HTO_classification)){
        
        sobj.subset <- subset(sobj, subset = HTO_classification == treatment)
        dim(sobj.subset)
        
        # calculate % mitocondrial counts
        sobj.subset[["percent.mt"]] <- PercentageFeatureSet(sobj.subset, pattern = "^MT-")
        
        # Visualize QC metrics as a violin plot. This can be used to decide values for filtering
        p <- VlnPlot(sobj.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        
        # save to a file
        png(paste(treatment, "QC_VlnPlot_Before_Filtering.png", sep = "_"),
            type = 'cairo',
            height = 4,
            width = 5,
            units = 'in',
            res = 200
        )
        print(p)
        dev.off()
        
        # get the counts of the largest expressed genes
        sobj.subset$largest_count <- apply(sobj.subset@assays$RNA@counts, 2, max)
        sobj.subset$largest_index <- apply(sobj.subset@assays$RNA@counts, 2, which.max)
        sobj.subset$largest_gene <- rownames(sobj.subset)[sobj.subset$largest_index]
        sobj.subset$percent.Largest.Gene <- sobj.subset$largest_count/sobj.subset$nCount_RNA * 100
        
        VlnPlot(sobj.subset, features=c("nCount_RNA","percent.mt", "percent.Largest.Gene")) #+ scale_y_log10()
        
        qc.metrics <- as_tibble(
            sobj.subset[[]],
            rownames="Cell.Barcode"
        ) 
        
        head(qc.metrics)
        
        
        # qc.metrics %>%
        #     arrange(percent.mt) %>%
        #     ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
        #     geom_point() + 
        #     scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
        #     ggtitle("Example of plotting QC metrics") +
        #     geom_hline(yintercept = 750) +
        #     geom_hline(yintercept = 2000) 
        
        # plot with log10 scale
        # p <- qc.metrics %>%
        #     arrange(percent.mt) %>%
        #     ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
        #     geom_point(size=0.7) + 
        #     scale_color_gradientn(colors=c("black","blue","green2","red","yellow")) +
        #     ggtitle(paste("QC metrics for", treatment)) +
        #     geom_hline(yintercept = 750) +
        #     geom_hline(yintercept = 2000) +
        #     scale_x_log10() + scale_y_log10()
        # 
        # png(paste(treatment, "Feature_QC_plot.png", sep = "_"),
        #     type = 'cairo',
        #     height = 5,
        #     width = 5,
        #     units = 'in',
        #     res = 300
        # )
        # print(p)
        # dev.off()
            
        # get largest gene list
        largest_gene_list <- qc.metrics %>%
            group_by(largest_gene) %>%
            count() %>%
            arrange(desc(n)) 
        
        largest_genes_to_plot <- largest_gene_list %>%
            filter(n>140) %>%
            pull(largest_gene)
        
        p <- qc.metrics %>%
            filter(largest_gene %in% largest_genes_to_plot) %>%
            mutate(largest_gene=factor(largest_gene, levels=largest_genes_to_plot)) %>%
            arrange(largest_gene) %>%
            ggplot(aes(x=log10(nCount_RNA), y=log10(nFeature_RNA), colour=largest_gene)) +
            geom_point(size=1) +
            scale_colour_manual(values=c("grey",RColorBrewer::brewer.pal(9,"Set1")))
        
        png(paste(treatment, "Largest_Gene_Scatterplot.png", sep = "_"),
            type = 'cairo',
            height = 5,
            width = 5,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        p <- qc.metrics %>%
            ggplot(aes(percent.Largest.Gene)) + 
            geom_histogram(binwidth = 0.7, fill="yellow", colour="black") +
            ggtitle("Distribution of Percentage Largest Gene") +
            geom_vline(xintercept = 10)
        
        png(paste(treatment, "Distribution_of_Largest_Gene_Histogram.png", sep = "_"),
            type = 'cairo',
            height = 5,
            width = 5,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        # based on the QC plots, we set the threshold for these parameters
        if (treatment == "MCF7RBposDMSO"){
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 7000
            nCount_RNA_min <- 4500
            nCount_RNA_max <- 31000
            percentMt <- 6
        } else if (treatment == "MCF7RBposPalbo"){
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 6500
            nCount_RNA_min <- 4500
            nCount_RNA_max <- 31000
            percentMt <- 7
        } else if (treatment == "MCF7RBneg"){
            nFeature_RNA_min <- 2200
            nFeature_RNA_max <- 7000
            nCount_RNA_min <- 4500
            nCount_RNA_max <- 31000
            percentMt <- 6
        } else if (treatment == "T47DRBposDMSO"){
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 7000
            nCount_RNA_min <- 2000
            nCount_RNA_max <- 31000
            percentMt <- 11
        } else if (treatment == "T47DRBposPalbo"){
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 7000
            nCount_RNA_min <- 4500
            nCount_RNA_max <- 31000
            percentMt <- 12
        } else if (treatment == "T47DRBneg"){
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 6000
            nCount_RNA_min <- 4500
            nCount_RNA_max <- 30000
            percentMt <- 10
        } else if (treatment == "HCC1806DMSO"){
            nFeature_RNA_min <- 3000
            nFeature_RNA_max <- 6000
            nCount_RNA_min <- 4500
            nCount_RNA_max <- 30000
            percentMt <- 6.0 
        } else if (treatment == "HCC1806Palbo"){
            nFeature_RNA_min <- 3000
            nFeature_RNA_max <- 6500
            nCount_RNA_min <- 10000
            nCount_RNA_max <- 32000
            percentMt <- 6.5
        } else if (treatment == " MiaPaCasgctrlDMSO"){
            nFeature_RNA_min <- 3000
            nFeature_RNA_max <- 6000
            nCount_RNA_min <- 10000
            nCount_RNA_max <- 30000
            percentMt <- 6
        } else if (treatment == "MiaPaCasgctrlPalbMRTX"){
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 5500
            nCount_RNA_min <- 4000
            nCount_RNA_max <- 25000
            percentMt <- 7
        } else if (treatment == "MiaPaCasgRBDMSO"){ # need to discuss with Erik for appropriate percentMt values
            nFeature_RNA_min <- 2500
            nFeature_RNA_max <- 6000
            nCount_RNA_min <- 3000
            nCount_RNA_max <- 30000
            percentMt <- 8
        } else if (treatment == "MiaPaCasgRBPalbMRTX"){ # need to discuss with Erik for appropriate percentMt values
            nFeature_RNA_min <- 2000
            nFeature_RNA_max <- 5500
            nCount_RNA_min <- 3000
            nCount_RNA_max <- 25000
            percentMt <- 8
        } else if (treatment == "MiaPaCasgRBPalbMRTX"){ # need to discuss with Erik for appropriate percentMt values
            nFeature_RNA_min <- 2500
            nFeature_RNA_max <- 6500
            nCount_RNA_max <- 
                nCount_RNA_min <- 
                percentMt <- 14
        } 
        
        # filter based on the parameters above to remove counts coming from died or lysed cells
        sobj.subset <- subset(sobj.subset, subset = nFeature_RNA > nFeature_RNA_min & nCount_RNA > nCount_RNA_min & 
                                  nCount_RNA < nCount_RNA_max & nFeature_RNA < nFeature_RNA_max &
                                  percent.mt < percentMt & percent.Largest.Gene < 10) 
        
        # Visualize QC metrics as a violin plot
        p <- VlnPlot(sobj.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        
        png(paste(treatment, "QC_VlnPlot_After_Filtering.png", sep = "_"),
            type = 'cairo',
            height = 4,
            width = 5,
            units = 'in',
            res = 200
        )
        print(p)
        dev.off()
        
        # perform normalization and scaling
        sobj.subset <- SCTransform(sobj.subset, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
        #sobj.subset <- NormalizeData(sobj.subset, normalization.method = "LogNormalize")
        #sobj.subset <- NormalizeData(sobj.subset, normalization.method = "CLR", margin = 2)
        
       #  ggplot(mapping = aes(sobj.subset@assays$RNA@data["GAPDH",])) + 
       #      geom_histogram(binwidth = 0.05, fill="yellow", colour="black") + 
       #      ggtitle("GAPDH expression")
       #  as.tibble(
       #      sobj.subset@assays$RNA@data[,1:100]
       #  ) %>%
       #      pivot_longer(
       #          cols=everything(),
       #          names_to="cell",
       #          values_to="expression"
       #      )  %>%
       #      ggplot(aes(x=expression, group=cell)) +
       #      geom_density() +
       #      coord_cartesian(ylim=c(0,0.6), xlim=c(0,3))
       #  
       #  tibble(
       #      pc95 = apply(sobj.subset[["RNA"]]@data,2,quantile,0.95),
       #      measured = apply(sobj.subset[["RNA"]]@data,2,function(x)(100*sum(x!=0))/length(x))
       #  ) -> normalisation.qc
       #  
       # # normalisation.qc %>% 
       #      ggplot(aes(x=measured,y=pc95))+
       #      geom_point()+
       #      ggtitle("Normalisation of data")
       #  
        # sobj.subset <- CellCycleScoring(sobj.subset, s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE)
        
        # as_tibble(sobj.subset[[]]) %>%
        #     ggplot(aes(Phase)) + geom_bar()
        # 
        # as_tibble(sobj.subset[[]]) %>%
        #     ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
        #     geom_point() +
        #     coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))
        # 
        
        # sobj.subset <- FindVariableFeatures(
        #     sobj.subset, 
        #     selection.method = "vst", 
        #     nfeatures=500
        # )
        # 
        # sobj.subset <- ScaleData(sobj.subset,features=rownames(sobj.subset))
        # sobj.subset <- RunPCA(sobj.subset,features=VariableFeatures(sobj.subset))
        # 
        #DimPlot(sobj.subset,reduction="pca")
        
        #DimHeatmap(sobj.subset,dims=1:15, cells=500)
        
        
        # perform linear dimensionality reduction
        sobj.subset <- RunPCA(sobj.subset, verbose = FALSE)
        sobj.subset <- RunUMAP(sobj.subset, dims = 1:30, verbose = FALSE)

        sobj.subset <- FindNeighbors(sobj.subset, dims = 1:30, verbose = FALSE)
        sobj.subset <- FindClusters(sobj.subset, verbose = FALSE)
        DimPlot(sobj.subset, label = TRUE) + NoLegend()

        VizDimLoadings(sobj.subset, dims = 1:2, reduction = "pca")
        
        # trajectory detection
        expression <- t(as.matrix(sobj.subset@assays$RNA@data))
        
        group_name <- sobj.subset@meta.data$HTO_classification 
        group_name <- as.factor(group_name)
        
        space <- reduce_dimensionality(expression, dist = "spearman", ndim = 3)
        
        traj <- infer_trajectory(space)
        
        p <- draw_trajectory_plot(
            space, 
            progression_group = group_name,
            path = traj$path,
            contour = TRUE
        )
        
        png(paste(treatment, "TrajectoryPlot.png", sep = "_"),
            type = 'cairo',
            height = 4,
            width = 5,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        
        # trajectory heatmap for select genes
        genes.of.interest <- c("CDK1", "PCNA", "CCND1", "CCNB1", "MCM2",  "RB1", "CCNE1", "CDT1", 
                               "RBL1", "RBL2", "E2F3", "E2F6", "CDKN1A", "CDKN1B")
        
        e2fTargetGenes <- c("AURKA", "BRCA2", "CCP110", "CENPE", "CKS2", "DCLRE1B", "DNMT1", "DONSON", "EED", "GINS1", 
                            "GINS4", "H2AFZ", "LIG1", "MAD2L1", "MCM2", "MCM4", "MCM5", "MCM7", "MELK", "MMS22L", 
                            "NAA38", "NASP", "NUDT21", "NUP205", "ORC6", "PCNA", "PLK4", "POLE", "PRIM2", "RAD51AP1", 
                            "RFC2", "RPA2", "RPA3", "SUV39H1", "TMPO", "UBE2T", "WDR90", "CDK1", "MCM3", "TOP2A", 
                            "MCM6", "BIRC5", "CCNB2", "RRM2", "HMGB2", "BUB1B", "RFC3", "EZH2", "CHEK1", "SMC4", 
                            "MKI67", "CDC20", "PLK1", "KIF2C", "DLGAP5", "AURKB", "CDC25A", "TRIP13", "H2AFX", "HMMR", 
                            "E2F8", "BRCA1", "MYBL2", "POLD1", "RACGAP1", "CKS1B", "KPNA2", "MSH2", "CDKN3", "ATAD2", 
                            "RPA1", "STMN1", "TIPIN", "TK1", "CDCA8", "ESPL1", "NCAPD2", "RANBP1", "MRE11", "KIF4A", 
                            "LMNB1", "KIF22", "UNG", "SMC1A", "CCNE1", "CDCA3", "ASF1B", "POLA2", "TIMELESS", 
                            "HELLS", "UBE2S", "PRKDC", "RAN", "USP1", "SPAG5", "POLD3", "DUT", "TACC3", "KIF18B", 
                            "CDC25B", "SRSF1", "GINS3", "NOLC1", "SLBP", "CHEK2", "SPC25", "BARD1", "DCTPP1", "SMC3", 
                            "RNASEH2A", "DEK", "CENPM", "RAD51C", "CBX5", "RFC1", "POLD2", "DSCC1", "ILF3", "DEPDC1", 
                            "DCK", "CDKN2C", "MYC", "TCF19", "RAD1", "LBR", "NBN", "PTTG1", "UBR7", "POLE4", "TUBG1",
                            "CTCF", "CNOT9", "TUBB", "SMC6", "ZW10", "PA2G4", "SSRP1", "NAP1L1", "ANP32E", "HMGB3", "IPO7", 
                            "RAD21", "CDK4", "CDKN1A", "BRMS1L", "CTPS1", "RAD50", "TRA2B", "CSE1L", "PAICS", "STAG1", 
                            "LUC7L3", "PPM1D", "NME1", "SRSF2", "XPO1", "HNRNPD", "PMS2", "ASF1A", "EXOSC8", "MLH1", 
                            "NUP107", "ORC2", "TP53", "TFRC", "HMGA1", "PSIP1", "DDX39A", "SNRPB", "CDKN1B", "MTHFD2", 
                            "WEE1", "PRDX4", "PHF5A", "TBRG4", "SHMT1", "PRPS1", "DIAPH3", "NUP153", "PSMC3IP",  "XRCC6", 
                            "PNN", "HUS1", "RBBP7", "PDS5B", "NOP56", "MXD3", "PPP1R8", "GSPT1",  "AK2", 
                            "CIT", "ING3", "JPT1", "POP7", "SYNCRIP", "EIF2S1", "LYAR", "PAN2", "SPC24") # "CDKN2A",
        
        expr_sel1 <- expression[, genes.of.interest]
        
        #e2fTargetGenes[!e2fTargetGenes %in% colnames(expression)]
        
        expr_sel <- expression[, e2fTargetGenes]
        
        h <- ncol(expr_sel) * 0.14 + 2
        w <- nrow(expr_sel) * 0.015 + 1
        
        #modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
        #p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = TRUE)
        p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, show_labels_row = TRUE)
        png(paste(treatment, "Trajectory_Heatmap_for_E2F_Target_Genes.png", sep = "_"),
            type = 'cairo',
            height = h,
            width = w,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        #modules <- extract_modules(scale_quantile(expr_sel1), traj$time, verbose = FALSE)
        #p <- draw_trajectory_heatmap(expr_sel1, traj$time, group_name, modules, show_labels_row = TRUE)
        p <- draw_trajectory_heatmap(expr_sel1, traj$time, group_name, show_labels_row = TRUE)
        h <- ncol(expr_sel1) * 0.14 + 2
        w <- nrow(expr_sel1) * 0.015 + 1
        png(paste(treatment, "Trajectory_Heatmap_for_Select_Genes.png", sep = "_"),
            type = 'cairo',
            height = h,
            width = w,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        # Make a line plot to show changes of scaled gene expression level. 
        # Input data were smoothed multiple times by repeatedly averaging values over sliding windows
        x <- expr_sel[order(traj$time), ]
        x <- scale_quantile(x, outlier_cutoff = 0.05)
        x.smoothed <- apply(x, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed$time <- 1:nrow(x.smoothed)
        
        dd = melt(x.smoothed, id=c("time"))
        
        colnames(dd) <- c('Time', 'Gene', 'Relative Expression')
        
        
        p2 <- ggplot(dd) + geom_line(aes(x=Time, y=`Relative Expression`, colour=Gene)) +
            #scale_colour_manual(values=c("red","green","blue", "purple")) +
            labs(title=treatment) +
            theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
        
        png(paste(treatment, "Trajectory_Line_Plot_For_E2F_Target_Genes.png", sep = "_"),
            type = 'cairo',
            height = 3, #h,
            width = 10,
            units = 'in',
            res = 200
        )
        print(p2)
        dev.off()
        
        # for genes of interest
        x <- expr_sel1[order(traj$time), ]
        x <- scale_quantile(x, outlier_cutoff = 0.05)
        x.smoothed <- apply(x, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed$time <- 1:nrow(x.smoothed)
        
        dd = melt(x.smoothed, id=c("time"))
        
        colnames(dd) <- c('Time', 'Gene', 'Relative Expression')
        p2 <- ggplot(dd) + geom_line(aes(x=Time, y=`Relative Expression`, colour=Gene)) +
            #scale_colour_manual(values=c("red","green","blue", "purple")) +
            labs(title=treatment) +
            theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8)) #+ NoLegend()
        
        png(paste(treatment, "Trajectory_Line_Plot_For_Select_Genes.png", sep = "_"),
            type = 'cairo',
            height = 5, #h,
            width = 10,
            units = 'in',
            res = 200
        )
        print(p2)
        dev.off()
        
        # # retrieve genes within each module
        # module_ids <- unique(modules$module)
        # for (mid in module_ids){
        #     features <- modules$feature[modules$module == mid]
        #     
        #     write.table(x = features,
        #                 file = paste(treatment, '_Genes_In_Module_', mid, '.txt', sep = ''),
        #                 quote = FALSE,
        #                 col.names = FALSE,
        #                 row.names = FALSE
        #     )
        # }
        
        # identify genes with expression inversely correlated with E2F target genes
        
        sobj.subset <- FindVariableFeatures(sobj.subset, selection.method = "vst", nfeatures = 2000)
        var.genes <- sobj.subset@assays$SCT@var.features
        var.data <- as.matrix(sobj.subset@assays$RNA@data[var.genes, ])
        
        #scaledData <- as.data.frame(sobj.subset@assays$SCT@scale.data)
        #corMatrix <- cor(t(scaledData), use="c")
        corMatrix <- cor(t(var.data), use="c")
        
        # get a subset of corMatrix to include only the E2F target genes
        corMatrixE2fTargets <- corMatrix[rownames(corMatrix) %in% e2fTargetGenes, colnames(corMatrix) %in% e2fTargetGenes]
        
        p <- pheatmap(corMatrixE2fTargets, fontsize_row = 8, fontsize_col = 8)
        
        png(paste(treatment, "Cor_Heatmap_for_E2F_Target_Genes.png", sep = "_"),
            type = 'cairo',
            height = 12,
            width = 12,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        
        # based on this heatmap, MKI67 forms a cluster with a few highly correlated genes. We will use this gene to get reversely 
        # correlated genes
        
        mki67CorGenes <- corMatrix["MKI67", ]
        mki67RevCorGenes <- names(mki67CorGenes[mki67CorGenes < -0.1])
        mki67CorGenes <- names(mki67CorGenes[mki67CorGenes > 0.6])
        
        if (length(mki67RevCorGenes) == 0){
            print(paste("No reversely correlated gene found for", treatment))
            next
        }
        expr_sel2 <- expression[, c(mki67CorGenes, mki67RevCorGenes)]
        h <- ncol(expr_sel2) * 0.14 + 2
        w <- nrow(expr_sel2) * 0.015 + 1
        modules <- extract_modules(scale_quantile(expr_sel2), traj$time, verbose = FALSE)
        p <- draw_trajectory_heatmap(expr_sel2, traj$time, group_name, modules, show_labels_row = TRUE)
        
        png(paste(treatment, "Trajectory_Heatmap_for_Inverse_Correleated_E2F_Target_Genes_.png", sep = "_"),
            type = 'cairo',
            height = h,
            width = w,
            units = 'in',
            res = 300
        )
        print(p)
        dev.off()
        
        # for inverse correlated genes
        x <- expr_sel2[order(traj$time), ]
        x <- scale_quantile(x, outlier_cutoff = 0.05)
        x.smoothed <- apply(x, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
        x.smoothed <- as.data.frame(x.smoothed)
        x.smoothed$time <- 1:nrow(x.smoothed)
        
        dd = melt(x.smoothed, id=c("time"))
        
        colnames(dd) <- c('Time', 'Gene', 'Relative Expression')
        p2 <- ggplot(dd) + geom_line(aes(x=Time, y=`Relative Expression`, colour=Gene)) +
            labs(title=treatment) +
            theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8))
        
        png(paste(treatment, "Trajectory_Line_Plot_For_Inversely_Corr_Genes.png", sep = "_"),
            type = 'cairo',
            height = 5, #h,
            width = 10,
            units = 'in',
            res = 200
        )
        print(p2)
        dev.off()
        
        
        gc()
    }
    
}
