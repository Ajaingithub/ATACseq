## This function is downstream analysis after filtering the peaks and removing the low quality and blacklisted region
# Args:
# cnts_path: this is a filter count matrix
# sampleinfo: sample metadata
# saveDir: saving directory
# group_level: Sorting of the groups in the order like Young, Old, Older
# label: In PCA ggplot label column for each points
# color: In PCA ggplot color column for each points
# shape: In PCA ggplot shape column for each points

PCA_plot <- function(cnts_path, sampleinfo, saveDir, group_level, label, color="Group", shape = NULL){
  message("Loading the counts and the target file\n")
  cnts_forDESeq <- read.table(cnts_path, sep = "\t", header = T)
  targes_forDESeq <- sampleinfo
  targes_forDESeq$Group <- tolower(targes_forDESeq$Group)
  print(table(targes_forDESeq$Group))
  counts_sum <- as.data.frame(as.matrix(colSums(cnts_forDESeq)))
  
  # while loading the cnts the colname change from -  to ., so making changes in the bam location in the targes for Deseq file also
  targes_bamlocation <- gsub("-",".",basename(targes_forDESeq$SampleID2))
  targes_forDESeq <- targes_forDESeq[match(rownames(counts_sum), targes_bamlocation),]
  
  message("Checking if both the counts and target file are in order\n")
  condition <- all(rownames(counts_sum) == targes_forDESeq$SampleID2)
  stopifnot(condition == "TRUE")
  targes_forDESeq$label_Group <- paste(targes_forDESeq[,label], targes_forDESeq[,color],sep = ":")
  targes_forDESeq$counts <- counts_sum$V1
  targes_forDESeq$Group<- factor(targes_forDESeq$Group, levels=group_level)
  condition <- all(targes_bamlocation == colnames(cnts_forDESeq))
  stopifnot(condition == TRUE)
  
  message("Making Deseq format for the counts")
  dds<-DESeqDataSetFromMatrix(countData = cnts_forDESeq,
                              colData = targes_forDESeq,
                              design = formula(~ Group))
  
  message("Filtering the data\n")
  keep <- filterByExpr(dds, group=dds$Group)
  dds <- dds[keep,]
  print(table(keep))
  
  message("Saving the dds to be used in the differential\n")
  saveRDS(dds, paste(saveDir,"dds.RDS",sep = ""))
  
  message("Running rlog transformation\n")
  set.seed(2021)
  rld <- rlog(dds) # rlog() may take a long time with 50 or more samples, that is why using vst
  saveRDS(rld, paste(saveDir,"dds_rlog.RDS",sep = ""))
  #rld <- readRDS(paste(saveDir,"dds_rlog.RDS",sep = ""))
  
  message("Running PCA")
  PCA <- DESeq2::plotPCA(object = rld, 
                         intgroup = colnames(targes_forDESeq), 
                         returnData=TRUE, ntop=5000)
  
  percentVar <- round(100 * attr(PCA, "percentVar"))
  percentVar
  
  PC1 <- ggplot(PCA, aes_string("PC1", "PC2", label="label_Group", color=color, shape = shape)) + 
    geom_point(size = 3) +
    geom_text_repel(size=4) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle("PCA EBV") +
    theme(plot.title = element_text(hjust = 0.5), 
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.minor = element_line(colour = "grey"),
          panel.grid.major = element_line(colour = "grey"))
  
  message("Saving PCA plot\n")
  pdf(paste(saveDir,"PCA.pdf",sep = ""), width = 15, height = 10)
  print(PC1)
  dev.off()
  
  message("Running UMAP\n")
  set.seed(2021)
  rld_cnts <- as.data.frame(rld@assays@data[[1]])
  cnts_forDESeq_umap <- umap(t(rld_cnts)) # others are default
  cnts_forDESeq_umap_layout <- as.data.frame(cnts_forDESeq_umap$layout)
  colnames(cnts_forDESeq_umap_layout) <- c("UMAP1","UMAP2")
  condition <- all(rownames(cnts_forDESeq_umap_layout) == targes_forDESeq$SampleID2)
  stopifnot(condition == TRUE)
  
  cnts_forDESeq_umap_layout_2 <- cbind(cnts_forDESeq_umap_layout,targes_forDESeq)
  p <- ggplot(cnts_forDESeq_umap_layout_2,
              aes_string("UMAP1","UMAP2", label="label_Group", color = color, shape=shape)) + 
    geom_point(size=4) + 
    geom_text_repel(size=4) + 
    theme_bw()
  
  message("Saving UMAP plot\n")
  pdf(paste(saveDir,"UMAP.pdf",sep = ""), width = 15, height = 10)
  print(p)
  dev.off()
  
  message("Please find the PCA and UMAP at this location ",saveDir,"\n")
  
}




