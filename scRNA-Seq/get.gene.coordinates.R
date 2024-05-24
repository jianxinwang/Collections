library(dplyr)
up.peaks.519 <- read.csv("/mnt/Jason/ATAC-Seq/output/ATAC-Seq_RBdCDK_plusDox_vs_minusDox_DE_Peaks_Up_Peaks_519_Homer_annotated.txt", sep = "\t")

# Exclude peaks in intergenic regions
up.peaks.519 <- up.peaks.519[up.peaks.519$Annotation != 'Intergenic', ]

gtf <- read.csv("/mnt/Jason/util/hg38.gtf", sep = '\t', header = FALSE)

gtf <- gtf[!grepl('_', gtf$V1), ]

gtf$refSeq <- sub(".*?(NM_\\d+).*", "\\1", gtf$V9)

genes <- gtf %>% group_by(V1, V9) %>% summarise(start=min(V4), stop=max(V5), .groups = 'drop')
genes$refSeq <- sub(".*?(NM_\\d+).*", "\\1", genes$V9)

merged <- merge(genes, up.peaks.519, by.x = "refSeq", by.y = 'Nearest.Refseq')

merged$id <- paste0(merged$refSeq, '_', merged$V1, ':',  merged$start, '-', merged$stop)

res <- merged[, c("V1", "start", "stop", "id")]
res <- res[!duplicated(res),]


write.table(x = res,
            file = "ATAC_up_peaks_519_whole_genes.bed",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t")
