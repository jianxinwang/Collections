library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(Seurat)
library(patchwork)
library(SCORPIUS)
library(ggplot2)
library(reshape2)
library(sctransform)
library(stats)
library(ggforce)

date <- Sys.Date()

smoothByWindow <- function(x){
    res <- NULL
    y <- c(x[(length(x) - 50):length(x)], x, x[1:50])
    
    for (i in 1:(length(x))){
        new.i <- mean(y[i:(i + 100)])
        res <- c(res, new.i)
    }
    res
}

fitSS <- function(xy,
                  a0=mean(xy[,1]),
                  b0=mean(xy[,2]),
                  r0 = mean(sqrt((xy[,1]-a0)^2 + (xy[,2]-b0)^2)),
                  ...){
  SS <- function(abr){
    sum((abr[3] - sqrt((xy[,1]-abr[1])^2 + (xy[,2]-abr[2])^2))^2)
  }
  optim(c(a0,b0,r0), SS, ...)
}

circlexy <- function(xyr, n=180){
  theta = seq(0,2*pi,len=n)
  cbind(xyr[1] + xyr[3]*cos(theta),
        xyr[2] + xyr[3]*sin(theta)
  )
}

dir <- "/Volumes/PersonalizedMed$/Jason/"

metaDataMCF7_T47D <- read.csv("/Volumes/PersonalizedMed$/Jason/scRNA/data/RS-03616030_ERpos-GEX_RS-03616038_ERpos-FB/metadata.txt", sep = "\t", row.names = 1)
metaDataHcc1806_MiaPaCa <- read.csv("/Volumes/PersonalizedMed$/Jason/scRNA/data/RS-03616031_MiaHCC-GEX_RS-03616039_MiaHCC-FB/metadata.txt", sep = "\t", row.names = 1)

expDataMCF7_T47D <- Read10X(data.dir = "/Volumes/PersonalizedMed$/Jason/scRNA/data/RS-03616030_ERpos-GEX_RS-03616038_ERpos-FB/")
expDataHCC1806_MiaPaCa <- Read10X(data.dir = "/Volumes/PersonalizedMed$/Jason/scRNA/data/RS-03616031_MiaHCC-GEX_RS-03616039_MiaHCC-FB//")

sobj <- CreateSeuratObject(counts = expDataMCF7_T47D$`Gene Expression`, project = "CellCycle") 
sobj[["HTO"]] <- CreateAssayObject(counts = expDataMCF7_T47D$`Antibody Capture`) 
sobj <- NormalizeData(sobj, assay = "HTO", normalization.method = "CLR")

sobj <- HTODemux(sobj, assay = "HTO", positive.quantile = 0.99)

table(sobj$HTO_classification.global)

# Visualize enrichment for selected HTOs with ridge plots
Idents(sobj) <- "HTO_maxID"
RidgePlot(sobj, assay = "HTO", features = rownames(sobj[["HTO"]])[1:2], ncol = 2)

# Visualize pairs of HTO signals to confirm mutual exclusivity in singlets
FeatureScatter(sobj, feature1 = "B0252", feature2 = "B0256")

# Compare number of UMIs for singlets, doublets and negative cells
Idents(sobj) <- "HTO_classification.global"
VlnPlot(sobj, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

countList <- list()
countList[['MCF7_T47D']][['count']] <- expDataMCF7_T47D$`Gene Expression`
countList[['HCC1806_MiaPaCa']][['count']] <- expDataHCC1806_MiaPaCa$`Gene Expression`
countList[["MCF7_T47D"]][['metaData']] <- metaDataMCF7_T47D
countList[["HCC1806_MiaPaCa"]][['metaData']] <- metaDataHcc1806_MiaPaCa

setwd("/Volumes/PersonalizedMed$/Jason/scRNA/output/test")

sobjList <- list()

gc()

for (cellline in names(countList)){

    count <- countList[[cellline]][["count"]]
    metaData <- countList[[cellline]][["metaData"]]
    
    # ensure the correct ordering of rownames of metadata and column names of expression matrix
    metaData <- metaData[match(colnames(count), rownames(metaData)), ]
    all(colnames(count) == rownames(metaData))

    # Initialize the Seurat object with the raw (non-normalized data).
    sobj <- CreateSeuratObject(counts = count, project = "CellCycle", min.cells = 3, min.features = 200, meta.data = metaData)


    # calculate mitocondrial sequence content
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")

    treatments <- unique(metaData$CellType)
    treatments <- treatments[!treatments %in% c("Negative", "Multiplet")]
    
    for (treatment in treatments){
    
        sobj.subset <- subset(sobj, subset = CellType == treatment)
        
        # Visualize QC metrics as a violin plot. This can be used to decide values for filtering
        p <- VlnPlot(sobj.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        
        # clean treatment name
        treatment2 <- gsub('/', '_', treatment)
        treatment2 <- gsub(' ', '_', treatment2)
        
        # png(paste(treatment2, "QC_VlnPlot_Before_Filtering.png", sep = "_"),
        #     type = 'cairo',
        #     height = 4,
        #     width = 5,
        #     units = 'in',
        #     res = 200
        #     )
        # print(p)
        # dev.off()
        
        #plot1 <- FeatureScatter(sobj.subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
        #plot2 <- FeatureScatter(sobj.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
        #plot1 + plot2

   
        # filter to remove counts coming from died or lysed cells
        # sobj.subset <- subset(sobj.subset, subset = nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max & percent.mt < percentMt) 
        sobj.subset <- subset(sobj.subset, subset = nFeature_RNA > 2000 & nFeature_RNA < 7500 & percent.mt < 15) 
        
        
        # get final dataset size
        print(paste(treatment, dim(sobj.subset)))
        
        # Visualize QC metrics as a violin plot
        p <- VlnPlot(sobj.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        
        # png(paste(treatment2, "QC_VlnPlot_After_Filtering.png", sep = "_"),
        #     type = 'cairo',
        #     height = 4,
        #     width = 5,
        #     units = 'in',
        #     res = 200
        # )
        # print(p)
        # dev.off()
        
        # perform normalization and scaling
        sobj.subset <- SCTransform(sobj.subset, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

        
        # perform linear dimensionality reduction
        sobj.subset <- RunPCA(sobj.subset, verbose = FALSE)
        sobj.subset <- RunUMAP(sobj.subset, dims = 1:30, verbose = FALSE)
        
        sobj.subset <- FindNeighbors(sobj.subset, dims = 1:30, verbose = FALSE)
        sobj.subset <- FindClusters(sobj.subset, verbose = FALSE)
        p <- DimPlot(sobj.subset, label = TRUE) + NoLegend()
        
        outfile <- paste(date, treatment2, "DimPlot.png", sep = '_')
        outfile <- gsub(' ', '_', outfile)
                      
        png(outfile,
            width = 5,
            height = 5,
            units = 'in',
            res = 300)
        print(p)
        dev.off()
        
        #
        # trajectory detection
        #
        expression <- t(as.matrix(sobj.subset@assays$RNA@data))

        
        #cell.cycle.genes <- c("CCNB1", "CCND1", "CCNA2", "CCNE2") #, "PCNA", "CTD1", "CDKN1A", "CDKN1B")#, "MKI67", "GMNN")
        
        # cell.cycle.genes <- c("CCND1", "CCNA1", "CCNA2", "CCNE1", "CTD1", "CDKN1A", "CDKN1B", "MKI67", "GMNN")
        
        cell.cycle.genes <- c("CCND1", "CCNA2", "CCNE1", "CTD1", "CDKN1A", "CDKN1B", "MKI67", "GMNN", "RB1")
        
        exp.cell.cycle <- expression[, colnames(expression) %in% cell.cycle.genes]
      
        total <- apply(exp.cell.cycle, 1, sum)
        exp.cell.cycle <- exp.cell.cycle[total > 0, ]
        
        corMatrix <- cor(t(exp.cell.cycle), use="everything", method = "pearson")
        corMatrix <- abs(corMatrix)
        
        data.for.ccd <- 0.5 - abs(corMatrix)/2
        
        scaledCmd <- cmdscale(data.for.ccd)
        
        scaledCmd <- as.data.frame(scaledCmd)
        colnames(scaledCmd) <- c("CMD1", "CMD2")
        
        f = fitSS(scaledCmd)
        
        exp.cell.cycle <- as.data.frame(exp.cell.cycle)
        
        
        for (gene in names(exp.cell.cycle)){
          
          scaledCmd$color <- ifelse(exp.cell.cycle[[gene]] > mean(exp.cell.cycle[[gene]]), 'orange', 'blue') # arbitrary value. 
          #colorDf <- as.data.frame(apply(exp.cell.cycle, 2, function(x){(x-min(x))/(max(x)-min(x))}))
          #Exp <- colorDf[[gene]]
          p <- ggplot(scaledCmd, aes(x = CMD1, y = CMD2, color = color)) + geom_point(alpha = 0.4, size = 0.3) + labs(x = "CMD dimension 1", y = "CMD dimension 2") + 
            geom_circle(aes(x0=f$par[1], y0=f$par[2], r=f$par[3]), linetype = 2,
                        inherit.aes=FALSE) + coord_fixed() +  ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
               
          #p <- p + scale_color_gradient(low="blue", high="orange")
          
          png(paste(date, treatment2, gene, "ccd_CMD_Scatterplot.png", sep = '_'),
              width = 4,
              height = 4,
              units = 'in',
              res = 300
          )
          print(p)
          dev.off()
        }
        
        # center the data on the origin of the fitted circle
        scaledCmd$x1 <- scaledCmd$CMD1 - f$par[1]
        scaledCmd$y1 <- scaledCmd$CMD2 - f$par[2]
        
        scaledCmd$angle <- atan2(scaledCmd$y1, scaledCmd$x1)
        
        scaledCmd$angle2 <-  ifelse(scaledCmd$angle >= 0, scaledCmd$angle, 3.14*2 + scaledCmd$angle)
        
        
        scaledCmd <- scaledCmd[order(scaledCmd$angle2), ]
        
        data.for.heatmap <- exp.cell.cycle[match(rownames(scaledCmd), rownames(exp.cell.cycle)), ]
        data.for.heatmap <- apply(data.for.heatmap, 2, function(x){(x-min(x))/(max(x)-min(x))})
        
        p <- Heatmap(t(data.for.heatmap),
                     cluster_rows = F,
                     cluster_columns = F,
                     show_column_names = F,
                     name = 'Scaled\nExpression'
        )
  
  
        png(paste(date, treatment2, "ccd_CMD_heatmap.png", sep = '_'),
            width = 6,
            height = 2.5,
            units = 'in',
            res = 200
            )
        print(p)
        dev.off()
        
        exp.cell.cycle.smoothed <- apply(data.for.heatmap, 2, smoothByWindow)
        exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
        exp.cell.cycle.smoothed <- apply(exp.cell.cycle.smoothed, 2, smoothByWindow)
        exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
        exp.cell.cycle.smoothed <- apply(exp.cell.cycle.smoothed, 2, smoothByWindow)
        exp.cell.cycle.smoothed <- as.data.frame(exp.cell.cycle.smoothed)
        exp.cell.cycle.smoothed$time <- 1:nrow(exp.cell.cycle.smoothed)
 
        dd = melt(exp.cell.cycle.smoothed, id=c("time"))
       
        colnames(dd) <- c('Time', 'Gene', 'Relative Expression')
        p2 <- ggplot(dd) + geom_line(aes(x=Time, y=`Relative Expression`, colour=Gene)) + labs(title=treatment) +
        theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 6))
           
        png(paste(date, treatment2, "Trajectory_Line_Plot_For_Cell_Cycle_Genes.png", sep = "_"),
            type = 'cairo',
            height = 4,
            width = 6,
            units = 'in',
            res = 200
        )
        print(p2)
        dev.off()
        gc()
    }
    
}
    
    
    #      group_name <- sobj.subset@meta.data$CellType 
    #      group_name <- as.factor(group_name)
    #      
    #      space <- reduce_dimensionality(expression, dist = "spearman", ndim = 3)
    #      
    #      traj <- infer_trajectory(space)
    #      
    #      p <- draw_trajectory_plot(
    #          space, 
    #          progression_group = group_name,
    #          path = traj$path,
    #          contour = TRUE
    #      )
    #      
    #      png(paste(treatment2, "TrajectoryPlot.png", sep = "_"),
    #          type = 'cairo',
    #          height = 4,
    #          width = 5,
    #          units = 'in',
    #          res = 300
    #      )
    #      print(p)
    #      dev.off()
    #     
    #      gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
    #      gene_sel <- gimp[gimp$pvalue < 1e-05,]
    #      
    #      if (sum(!is.na(gene_sel$gene)) == 0){
    #          gene_sel <- gimp[1:50,]
    #      }
    #      
    #      expr_sel <- expression[,gene_sel$gene]
    #      
    #      h <- ncol(expr_sel) * 0.14 + 1
    #      w <- nrow(expr_sel) * 0.015 + 1
    #      
    #      modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
    #      p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = TRUE)
    #      p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, show_labels_row = TRUE)
    #      png(paste(treatment2, "Modulized_Trajectory_Heatmap.png", sep = "_"),
    #          type = 'cairo',
    #          height = h,
    #          width = w,
    #          units = 'in',
    #          res = 300
    #     )
    #      print(p)
    #      dev.off()
    #     
    #     
    #     # trajectory heatmap for select genes
    #     genes.of.interest <- c("CDK1", "PCNA", "CCND1", "CCNB1", "MCM2", "MKI67", "RB1", "CCNE1", "CDT1", 
    #                            "RBL1", "RBL2", "E2F3", "E2F1", "E2F6", "CDKN1A", "CDKN1B")
    #     
    #     e2fTargetGenes <- c("AURKA", "BRCA2", "CCP110", "CENPE", "CKS2", "DCLRE1B", "DNMT1", "DONSON", "EED", "GINS1", 
    #                         "GINS4", "H2AFZ", "LIG1", "MAD2L1", "MCM2", "MCM4", "MCM5", "MCM7", "MELK", "MMS22L", 
    #                         "NAA38", "NASP", "NUDT21", "NUP205", "ORC6", "PCNA", "PLK4", "POLE", "PRIM2", "RAD51AP1", 
    #                         "RFC2", "RPA2", "RPA3", "SUV39H1", "TMPO", "UBE2T", "WDR90", "CDK1", "MCM3", "TOP2A", 
    #                         "MCM6", "BIRC5", "CCNB2", "RRM2", "HMGB2", "BUB1B", "RFC3", "EZH2", "CHEK1", "SMC4", 
    #                         "MKI67", "CDC20", "PLK1", "KIF2C", "DLGAP5", "AURKB", "CDC25A", "TRIP13", "H2AFX", "HMMR", 
    #                         "E2F8", "BRCA1", "MYBL2", "POLD1", "RACGAP1", "CKS1B", "KPNA2", "MSH2", "CDKN3", "ATAD2", 
    #                         "RPA1", "STMN1", "TIPIN", "TK1", "CDCA8", "ESPL1", "NCAPD2", "RANBP1", "MRE11", "KIF4A", 
    #                         "LMNB1", "KIF22", "UNG", "SMC1A", "CCNE1", "CDCA3", "ASF1B", "POLA2", "TIMELESS", 
    #                         "HELLS", "UBE2S", "PRKDC", "RAN", "USP1", "SPAG5", "POLD3", "DUT", "TACC3", "KIF18B", 
    #                         "CDC25B", "SRSF1", "GINS3", "NOLC1", "SLBP", "CHEK2", "SPC25", "BARD1", "DCTPP1", "SMC3", 
    #                         "RNASEH2A", "DEK", "CENPM", "RAD51C", "CBX5", "RFC1", "POLD2", "DSCC1", "ILF3", "DEPDC1", 
    #                         "DCK", "CDKN2C", "MYC", "TCF19", "RAD1", "LBR", "NBN", "PTTG1", "UBR7", "POLE4", "TUBG1",
    #                         "CTCF", "CNOT9", "TUBB", "SMC6", "ZW10", "PA2G4", "SSRP1", "NAP1L1", "ANP32E", "HMGB3", "IPO7", 
    #                         "RAD21", "CDK4", "CDKN1A", "BRMS1L", "CTPS1", "RAD50", "TRA2B", "CSE1L", "PAICS", "STAG1", 
    #                         "LUC7L3", "PPM1D", "NME1", "SRSF2", "XPO1", "HNRNPD", "PMS2", "ASF1A", "EXOSC8", "MLH1", 
    #                         "NUP107", "ORC2", "TP53", "TFRC", "HMGA1", "PSIP1", "DDX39A", "SNRPB", "CDKN1B", "MTHFD2", 
    #                         "WEE1", "PRDX4", "PHF5A", "TBRG4", "SHMT1", "PRPS1", "DIAPH3", "NUP153", "PSMC3IP",  "XRCC6", 
    #                         "PNN", "HUS1", "RBBP7", "PDS5B", "NOP56", "MXD3", "PPP1R8", "GSPT1",  "AK2", 
    #                         "CIT", "ING3", "JPT1", "POP7", "SYNCRIP", "EIF2S1", "LYAR", "PAN2", "SPC24") # "CDKN2A",
    #     
    #     expr_sel1 <- expression[, genes.of.interest]
    #     
    #     e2fTargetGenes[!e2fTargetGenes %in% colnames(expression)]
    #     
    #      if (cellline == "HCC1806_MiaPaCa"){
    #          e2fTargetGenes <- e2fTargetGenes[!e2fTargetGenes %in% "CDKN2A"]
    #      }
    #     
    #     expr_sel <- expression[, e2fTargetGenes]
    #     
    #     h <- ncol(expr_sel) * 0.14 + 2
    #     w <- nrow(expr_sel) * 0.015 + 1
    #     modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
    #     p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = TRUE)
    #     
    #     png(paste(treatment2, "Trajectory_Heatmap_for_E2F_Target_Genes.png", sep = "_"),
    #         type = 'cairo',
    #         height = h,
    #         width = w,
    #         units = 'in',
    #         res = 300
    #     )
    #     print(p)
    #     dev.off()
    #     
    #     modules <- extract_modules(scale_quantile(expr_sel1), traj$time, verbose = FALSE)
    #     p <- draw_trajectory_heatmap(expr_sel1, traj$time, group_name, modules, show_labels_row = TRUE)
    #     
    #     h <- ncol(expr_sel1) * 0.14 + 2
    #     w <- nrow(expr_sel1) * 0.015 + 1
    #     png(paste(treatment2, "Trajectory_Heatmap_for_Select_Genes.png", sep = "_"),
    #         type = 'cairo',
    #         height = h,
    #         width = w,
    #         units = 'in',
    #         res = 300
    #     )
    #     print(p)
    #     dev.off()
    #     # Make a line plot to show changes of scaled gene expression level. 
    #     # Input data were smoothed multiple times by repeatedly averaging values over sliding windows
    #     x <- expr_sel[order(traj$time), ]
    #     x <- scale_quantile(x, outlier_cutoff = 0.05)
    #     x.smoothed <- apply(x, 2, smoothByWindow)
    #     x.smoothed <- as.data.frame(x.smoothed)
    #     x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
    #     x.smoothed <- as.data.frame(x.smoothed)
    #     x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
    #     x.smoothed <- as.data.frame(x.smoothed)
    #     x.smoothed$time <- 1:nrow(x.smoothed)
    #     
    #     dd = melt(x.smoothed, id=c("time"))
    # 
    #     colnames(dd) <- c('Time', 'Gene', 'Relative Expression')
    #     p2 <- ggplot(dd) + geom_line(aes(x=Time, y=`Relative Expression`, colour=Gene)) +
    #         #scale_colour_manual(values=c("red","green","blue", "purple")) +
    #         labs(title=treatment) +
    #         theme(plot.title = element_text(hjust = 0.5)) + NoLegend()
    #     
    #     png(paste(treatment2, "Trajectory_Line_Plot_For_E2F_Target_Genes.png", sep = "_"),
    #         type = 'cairo',
    #         height = 3, #h,
    #         width = 8,
    #         units = 'in',
    #         res = 300
    #     )
    #     print(p2)
    #     dev.off()
    #     
    #     # for genes of interest
    #     x <- expr_sel1[order(traj$time), ]
    #     x <- scale_quantile(x, outlier_cutoff = 0.05)
    #     x.smoothed <- apply(x, 2, smoothByWindow)
    #     x.smoothed <- as.data.frame(x.smoothed)
    #     x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
    #     x.smoothed <- as.data.frame(x.smoothed)
    #     x.smoothed <- apply(x.smoothed, 2, smoothByWindow)
    #     x.smoothed <- as.data.frame(x.smoothed)
    #     x.smoothed$time <- 1:nrow(x.smoothed)
    #     
    #     dd = melt(x.smoothed, id=c("time"))
    #     
    #     colnames(dd) <- c('Time', 'Gene', 'Relative Expression')
    #     p2 <- ggplot(dd) + geom_line(aes(x=Time, y=`Relative Expression`, colour=Gene)) +
    #         #scale_colour_manual(values=c("red","green","blue", "purple")) +
    #         labs(title=treatment) +
    #         theme(plot.title = element_text(hjust = 0.5), legend.text = element_text(size = 8)) #+ NoLegend()
    #     
    #     png(paste(treatment2, "Trajectory_Line_Plot_For_Select_Genes.png", sep = "_"),
    #         type = 'cairo',
    #         height = 3, #h,
    #         width = 8,
    #         units = 'in',
    #         res = 300
    #     )
    #     print(p2)
    #     dev.off()
    #     # retrieve genes within each module
    #     module_ids <- unique(modules$module)
    #     for (mid in module_ids){
    #         features <- modules$feature[modules$module == mid]
    #         
    #         write.table(x = features,
    #                     file = paste(treatment2, '_Genes_In_Module_', mid, '.txt', sep = ''),
    #                     quote = FALSE,
    #                     col.names = FALSE,
    #                     row.names = FALSE
    #                     )
    #     }
    #     
    #     # identify genes with expression inversely correlated with E2F target genes
    #     scaledData <- as.data.frame(sobj.subset@assays$SCT@scale.data)
    #     
    #     corMatrix <- cor(t(scaledData), use="c")
    #     
    #     # get a subset of corMatrix to include only the E2F target genes
    #     corMatrixE2fTargets <- corMatrix[rownames(corMatrix) %in% e2fTargetGenes, colnames(corMatrix) %in% e2fTargetGenes]
    #     
    #     p <- pheatmap(corMatrixE2fTargets)
    # 
    #     png(paste(treatment2, "Cor_Heatmap_for_E2F_Target_Genes.png", sep = "_"),
    #         type = 'cairo',
    #         height = 12,
    #         width = 12,
    #         units = 'in',
    #         res = 300
    #     )
    #     print(p)
    #     dev.off()
    #     
    #     
    #     # based on this heatmap, MKI67 forms a cluster with a few highly correlated genes. We will use this gene to get reversely 
    #     # correlated genes
    #      
    #     mki67CorGenes <- corMatrix["MKI67", ]
    #     mki67RevCorGenes <- names(mki67CorGenes[mki67CorGenes < -0.3])
    #     mki67CorGenes <- names(mki67CorGenes[mki67CorGenes > 0.6])
    #         
    #     if (length(mki67RevCorGenes) == 0){
    #         print(paste("No reversely correlated gene found for", treatment2))
    #         next
    #     }
    #     expr_sel <- expression[, c(mki67CorGenes, mki67RevCorGenes)]
    #     h <- ncol(expr_sel) * 0.14 + 2
    #     w <- nrow(expr_sel) * 0.015 + 1
    #     modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
    #     p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = TRUE)
    #     
    #     png(paste(treatment2, "Trajectory_Heatmap_for_Inverse_Correleated_E2F_Target_Genes_.png", sep = "_"),
    #         type = 'cairo',
    #         height = h,
    #         width = w,
    #         units = 'in',
    #         res = 300
    #     )
    #     print(p)
    #     dev.off()
    #     
    #     FeatureScatter(sobj.subset, feature1 = "rna_MKI67", feature2 = "rna_MKI67")
    #     
    # }
    # 
    # # combined treatment analysis
    # if (cellline == "HCC1806_MiaPaCa"){
    #     sobj.hcc1806.dmso.palbo <- subset(sobj, subset = CellType == "HCC1806 DMSO" |  CellType == "HCC1806 Palbo")
    #     sobj.MiaPaCa.dmso.palbo <- subset(sobj, subset = CellType == "MiaPaCa sg ctrl DMSO" | CellType == "MiaPaCa sg ctrl Palb/MRTX")
    #     
    #     sobj.hcc1806.dmso.palbo <- subset(sobj.hcc1806.dmso.palbo, subset = nFeature_RNA > 3000 & nFeature_RNA < 7000 & percent.mt < 7.5)
    #     sobj.MiaPaCa.dmso.palbo <- subset(sobj.MiaPaCa.dmso.palbo, subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 7)
    #     
    #     sobjList[["HCC1806"]] <- sobj.hcc1806.dmso.palbo
    #     sobjList[["MiaPaCa"]] <- sobj.MiaPaCa.dmso.palbo
    #              
    # } else if (cellline == "MCF7_T47D"){
    #     sobj.mcf7.dmso.palbo <- subset(sobj, subset = CellType == "MCF7 RB+ DMSO" | CellType == "MCF7 RB+ Palbo")
    #     sobj.t47d.dmso.palbo <- subset(sobj, subset = CellType == "T47D RB+ DMSO" | CellType == "T47D RB+ Palbo")
    #     
    #     sobj.mcf7.dmso.palbo <- subset(sobj.mcf7.dmso.palbo, subset = CellType == "MCF7 RB+ DMSO" | CellType == "MCF7 RB+ Palbo" | CellType == "MCF7 RB-")
    #     sobj.t47d.dmso.palbo <- subset(sobj.t47d.dmso.palbo, subset = CellType == "T47D RB+ DMSO" | CellType == "T47D RB+ Palbo" | CellType == "T47D RB-")
    #     
    #     sobjList[["MCF7"]] <- sobj.mcf7.dmso.palbo
    #     sobjList[["T47D"]] <- sobj.t47d.dmso.palbo
    # }
#    gc()
#}

# # combined dataset analysis
# for (cellline in names(sobjList)){
#     sobj.subset <- sobjList[[cellline]]
#     
#     # perform normalization and scalling
#     sobj.subset <- SCTransform(sobj.subset, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
#     
#     
#     # perform linear dimensionality reduction
#     sobj.subset <- RunPCA(sobj.subset, verbose = FALSE)
#     sobj.subset <- RunUMAP(sobj.subset, dims = 1:30, verbose = FALSE)
#     
#     sobj.subset <- FindNeighbors(sobj.subset, dims = 1:30, verbose = FALSE)
#     sobj.subset <- FindClusters(sobj.subset, verbose = FALSE)
#     DimPlot(sobj.subset, label = TRUE) + NoLegend()
#     
#     VizDimLoadings(sobj.subset, dims = 1:2, reduction = "pca")
#     
#     # trajectory detection
#     expression <- t(as.matrix(sobj.subset@assays$RNA@data))
#     group_name <- sobj.subset@meta.data$CellType 
#     group_name <- as.factor(group_name)
#     
#     space <- reduce_dimensionality(expression, dist = "spearman", ndim = 3)
#     
#     traj <- infer_trajectory(space)
#     
#     p <- draw_trajectory_plot(
#         space, 
#         progression_group = group_name,
#         path = traj$path,
#         contour = TRUE
#     )
#     
#     png(paste(cellline, "MultipleCellTypeTrajectoryPlot.png", sep = "_"),
#         type = 'cairo',
#         height = 4,
#         width = 5,
#         units = 'in',
#         res = 300
#     )
#     print(p)
#     dev.off()
#     
#     gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 8)
#     gene_sel <- gimp[gimp$pvalue < 1e-05,]
#     
#     if (sum(!is.na(gene_sel$gene)) == 0){
#         gene_sel <- gimp[1:50,]
#     }
#     
#     expr_sel <- expression[,gene_sel$gene]
#     
#     h <- ncol(expr_sel) * 0.14 + 1
#     w <- nrow(expr_sel) * 0.015 + 1
#     
#     modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
#     p <- draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = TRUE)
#     
#     png(paste(cellline, "Multiple_CellType_Modulized_Trajectory_Heatmap.png", sep = "_"),
#         type = 'cairo',
#         height = h,
#         width = w,
#         units = 'in',
#         res = 300
#     )
#     print(p)
#     dev.off()
#     
#     
#     
#}

