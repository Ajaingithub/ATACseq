########## Differential ##############
### Differentiation ####
# mainDir: The path where the 
# compare group should provide the detail of the group that needs to be compared
# main name is the name if the project for which session will be saved
# saveDir where all the plots and the files will be saved
# comparison the group that needs to be compared

ATAC_differential <- function(mainDir,saveDir, compare_group, main_name, comparison, pvalue = 0.05){
  `%notin%` <- Negate(`%in%`)
  message("Loading the Deseq file\n")
  dds <- readRDS(paste(mainDir,"dds.RDS",sep = ""))
  dds <- estimateSizeFactors(dds)
  
  jcp <- factor(dds[[compare_group]])
  design1 <- model.matrix(~ 0 + jcp)
  colnames(design1) <- levels(jcp)
  
  ## Comparing old_vzv_pre and young_vzv_pre
  message("\nMaking Contrast\n")
  
  # comparison[12],
  # comparison[13],comparison[14],comparison[15],comparison[16],comparison[17],comparison[18],
  # comparison[19],comparison[20],comparison[21],comparison[22],comparison[23],comparison[24],
  # comparison[25],comparison[26],comparison[27],comparison[28],comparison[29],comparison[30],
  # comparison[31],comparison[32],comparison[33],comparison[34],comparison[35],comparison[36],
  # comparison[37],comparison[38],comparison[39],comparison[40],comparison[41],comparison[42],
  # comparison[43],comparison[44],comparison[45],comparison[46],comparison[47],comparison[48],
  # comparison[49],comparison[50],
  cm <- makeContrasts(comparison[1],comparison[2],comparison[3],comparison[4],comparison[5],comparison[6],
                      comparison[7],comparison[8],comparison[9],comparison[10],comparison[11],  #comparison[12],
                      # comparison[13],comparison[14],comparison[15],comparison[16],comparison[17],comparison[18],
                      # comparison[19],comparison[20],comparison[21],#comparison[22],comparison[23],comparison[24],
                      # comparison[25],comparison[26],comparison[27],comparison[28],comparison[29],comparison[30],
                      # comparison[31],comparison[32],comparison[33],comparison[34],comparison[35],comparison[36],
                      # comparison[37],comparison[38],comparison[39],comparison[40],comparison[41],comparison[42],
                      # comparison[43],comparison[44],comparison[45],comparison[46],comparison[47],comparison[48],
                      # comparison[49],comparison[50],
                      levels = design1)
  cmdf<- as.data.frame(cm)[,colSums(is.na(as.data.frame(cm)))<nrow(as.data.frame(cm))]
  
  ## If you have not ran the ChromVar, before going for differential run this first
  peaks_filter <- as.data.frame(rownames(dds))
  peaks_filter_2 <- as.data.frame(gsub(":","\t",peaks_filter$`rownames(dds)`) %>% gsub("-","\t",.))
  colnames(peaks_filter_2) <- NULL
  
  write.table(peaks_filter_2,
              paste(saveDir,"peaks_filter_pre_2.txt",sep = ""),
              sep = "\t",
              quote = FALSE, 
              row.names = FALSE,
              col.names = FALSE)
  
  
  peaks <- getPeaks(paste(saveDir,"peaks_filter_pre_2.txt",sep = ""), sort_peaks=TRUE) # now we got the good quality peaks.
  
  
  message("\nPerforming GC Normalization")
  message("\nLoading the count matrix")
  cts <- counts(dds) # Getting to the counts
  FASTAfile <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
  FASTA <- FaFile(FASTAfile)
  open(FASTA)
  sequence <- getSeq(FASTA, peaks)
  GC_frequency <- letterFrequency(getSeq(FASTA, peaks), "GC")[,1]
  coordinate_Width <- width(peaks)
  peaks_df <- as.data.frame(peaks)
  coordinates <- paste(peaks_df$seqnames,":",peaks_df$start-1,"-",peaks_df$end,sep = "") 
  peaks_length_GC <- data.frame(coordinates, GC_frequency, coordinate_Width)
  peaks_length_GC$GC_fraction <- peaks_length_GC$GC_frequency/peaks_length_GC$coordinate_Width
  all(peaks_length_GC$coordinates == rownames(cts)) # getting all the match genes
  peaks_length_GC_2 <- peaks_length_GC[match(rownames(cts), peaks_length_GC$coordinates),]
  all(peaks_length_GC_2$coordinates == rownames(cts))
  fit <- cqn(cts, peaks_length_GC_2$GC_fraction, peaks_length_GC_2$coordinate_Width, sizeFactors = dds$sizeFactor)
  cqnOffset <-  fit$glm.offset
  cqnNormFactors <- exp(cqnOffset)
  #normalizationFactors(dds) <- cqnNormFactors
  
  message("\n Running Differential")
  dir.create(paste(saveDir,"Differential",sep = ""), showWarnings = FALSE)
  desl <- as.DGEList(dds)
  desl <- calcNormFactors(desl) # Normalizing for the cqn with size factor
  desl$offset <-  fit$glm.offset
  # Alternatively, gene-specific correction factors can be entered into the glm functions of edgeR as offsets. In the latter case, the offset matrix
  # will be assumed to account for all normalization issues, including sequencing depth and RNA composition.
  
  # Note that normalization in edgeR is model-based, and the original read counts are not themselves transformed. This means that users should not transform the read counts in any way
  #before inputing them to edgeR. For example, users should not enter RPKM or FPKM values to edgeR in place of read counts. 
  
  desl <- estimateGLMCommonDisp(desl, design = design1) # The pseudo-counts represent the equivalent counts would have been 
  # observed had the library sizes all been equal, assuming the fitted model. They are intended mainly for internal use in the edgeR pipeline.
  desl_Wnorm <- voomWithQualityWeights(desl,
                                       design1,
                                       normalize.method = 'quantile',
                                       plot = FALSE)
  
  desl_Wnorm_fit <- lmFit(desl_Wnorm,
                          design1,
                          method = 'robust',
                          maxit = 10000)
  
  message("saving the RDS lmfitted model\n")
  saveRDS(desl_Wnorm_fit,paste(saveDir,"Differential/desl_Wnorm_fit.RDS",sep = ""))
  
  desl_Wnorm_fit_2 <- contrasts.fit(desl_Wnorm_fit, cm)
  desl_Wnorm_fit_2 <- eBayes(desl_Wnorm_fit_2,  proportion = 1/10)
  results <- decideTests(desl_Wnorm_fit_2, p.value = 0.05)
  print(summary(results))
  
  dir.create(paste(saveDir,"Differential",sep = ""), showWarnings = FALSE)
  message("saving the fitted value result for\n")
  write.fit(fit = desl_Wnorm_fit_2, results = results,
            file = paste(saveDir,"Differential/All_region_wnorm_fit_",main_name,".txt",sep = ""),
            adjust="BH",sep = "\t")
  
  results.summary <- summary(results)
  
  message("saving the Differential gene number\n")
  write.table(t(results.summary[c(1, 3),]),
              file = paste(saveDir,"Differential/DE_genes_number.txt", sep = ""),
              quote = FALSE,
              row.names = TRUE,
              col.names = TRUE,
              sep = "\t")
  
  result_reshape <- reshape2::melt(t(results.summary[c(1, 3),]))
  colnames(result_reshape) <- c("comparison","Differential","Number_of_Genes")
  gsub(".*=","",result_reshape$comparison)
  h <- ggplot(result_reshape) + geom_bar(aes(x=comparison, y=Number_of_Genes, fill = Differential, color=Differential),
                                         stat = "identity", position = position_dodge()) +
    theme(plot.title = element_text(hjust = 0.5, size = 3),
          text = element_text(size = 13),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.minor = element_line(colour = "grey"),
          panel.grid.major = element_line(colour = "grey"))
  
  pdf(paste(saveDir,"Differential/Differential_gene_number.pdf",sep = ""), width = 8, height = 6)
  print(h)
  dev.off()
  
  message("Saving the MA plot")
  YF_MA_plot<- function(fit,file,px,py){
    xloc <- px   #position of the text for x axis
    yloc <- py   #position of the text for y axis
    pdf(paste(saveDir,"/Differential/",file,"_MA_plot.pdf", sep = ""),width = 10,height = 8)
    tops<-topTable(fit,
                   coef = file, 
                   number = Inf, 
                   sort.by = "none") #Extract a table of the top-ranked genes from a linear model fit
    tops_sig<-subset(tops,adj.P.Val<0.05)
    tops_sig_up<-subset(tops_sig,logFC>0)
    tops_sig_down<-subset(tops_sig,logFC < 0)
    tops_sig_mod<-tops_sig_up
    tops_sig_mod<-rbind(tops_sig_up,tops_sig_down)
    status <- row.names(tops)%in%row.names(tops_sig_mod)
    values <- c("FALSE","TRUE")
    col <- c("black","red")
    attr(status,"values") <- values
    attr(status,"col") <- col
    limma::plotMA(fit, coef = file, 
                  status = status,
                  cex=0.8,main = file,
                  legend=FALSE)
    message("Up=",nrow(tops_sig_up))
    message("Down=",nrow(tops_sig_down))
    text(xloc,yloc,paste0("Up=",nrow(tops_sig_up)))
    text(xloc,-yloc,paste0("Down=",nrow(tops_sig_down)))
    dev.off()
  }
  
  for(i in 1:length(colnames(cmdf))){
    YF_MA_plot(desl_Wnorm_fit_2, colnames(cmdf)[i], 8, 3)
  }
  
  message("saving the top changed genes")
  for(i in 1:length(colnames(cmdf))){
    tops<-topTable(desl_Wnorm_fit_2,
                   coef = colnames(cmdf)[i],
                   number = Inf,
                   sort.by = "none") #Extract a table of the top-ranked genes from a linear model fit
    all_region <- rownames(tops)
    tops_sig<-subset(tops,adj.P.Val<0.05)
    tops_sig_up<-subset(tops_sig,logFC>0)
    write.table(tops_sig_up, 
                paste(saveDir,"Differential/",colnames(cmdf)[i],"sig_up_foreground_all.txt",sep = ""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(rownames(tops_sig_up), 
                paste(saveDir,"Differential/",colnames(cmdf)[i],"sig_up_foreground.txt",sep = ""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    tops_sig_down<-subset(tops_sig,logFC < 0)
    write.table(tops_sig_down, 
                paste(saveDir,"Differential/",colnames(cmdf)[i],"sig_down_foreground_all.txt",sep = ""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(rownames(tops_sig_down), 
                paste(saveDir,"Differential/",colnames(cmdf)[i],"sig_down_foreground.txt",sep = ""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    write.table(rownames(tops_sig),
                paste(saveDir,"Differential/",colnames(cmdf)[i],"significant_peaks.txt",sep = ""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    all_region_except_differential <- all_region[which(all_region %notin% rownames(tops_sig))]
    write.table(all_region_except_differential,
                paste(saveDir,"Differential/",colnames(cmdf)[i],"non_significant_peaks_background.txt",sep = ""),
                quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  message("saving the session")
  save.session(paste(saveDir,"Differential/",main_name,".Rda",sep = ""))
  
  message(paste("Please go to this directory ",saveDir," to check the MA plots, DE gene number, top Genes, and bar plots", sep = ""))
}
