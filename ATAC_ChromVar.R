######## ChromVar ##########
# mainDir is the dir where the dds.RDS file has been saved
# saveDir is the dir where you want to save the file generated in the ChromVar
# sampleinfo is the metadata that is required 
# SampleID are the column name of the ID that are matched with the counts table header and the sampleinfo
# SampleDepth are the column name in the sampleinfo of the depth of the each sample that we get when we run the After Splot.
# SampleID_2 is the column in the sampleinfo you wants to put in the graph, this could be short name and easily identifiable
# col_arrange it is for the heatmap, so how do you want to arrange the group in the heatmap. Please provide the name like 
# col_arrange <- c("naïve","cm","em","temra","lat2","lyt1","lyt2"). Please copy paste the unique values from the sampleinfo

ATAC_ChromVar <- function(mainDir, saveDir, sampleinfo, SampleID, SampleDepth, SampleID_2, 
                          core=32, min_peaks=0.15, min_depth=1500, Group = "Group", col_arrange){
  message("putting multicore processor\n")
  register(MulticoreParam(core, progressbar = TRUE))
  
  message("Getting the peaks which are in FilterByExpr \n")
  dds <- readRDS(paste(mainDir,"dds.RDS",sep = ""))
  peaks_filter <- as.data.frame(rownames(dds))
  colnames(peaks_filter) <- NULL
  write.table(peaks_filter, paste(saveDir,"peaks_filter_pre.txt",sep = ""), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  message("Converting peaks to the bed files \n")
  setwd(saveDir)
  system ("sed \x27s/:/\t/g\x27 peaks_filter_pre.txt | sed \x27s/-/\t/g\x27 > peaks_filter_pre_2.txt")
  
  peaks <- getPeaks("peaks_filter_pre_2.txt", sort_peaks=TRUE) 
  # now we got the good quality peaks.
  # Ignore Warning we donot need the peaks of equal size
  peak_path <- paste(saveDir,"peaks_filter_pre.txt",sep = "")
  
  # we have saved the cnts data in the After_Splot.R as saveDir which will be considered here as mainDir
  counts_path <- paste(mainDir,"cnts_noblacklist_data.txt",sep = "")
  counts_filt_path <- paste(saveDir,"cnts_noblacklist_filtered_pre.txt",sep = "") # filter file will be saved in the saveDir
  system(paste("grep -wFf",peak_path,counts_path,">",counts_filt_path, sep = " "))
  system(paste("head -1",counts_path,"> header_cnts.txt",sep = " "))
  system("cat header_cnts.txt cnts_noblacklist_filtered_pre.txt > cnts_noblacklist_filtered_pre_2.txt")
  
  counts = read.table(paste(saveDir,"cnts_noblacklist_filtered_pre_2.txt",sep = ""), header = TRUE)
  ###########
  condition = all(colnames(counts) == sampleinfo[,SampleID])
  stopifnot(condition == "TRUE")
  my_counts = counts[,sampleinfo[,SampleID]]
  my_counts_matrix <- as.matrix(my_counts)
  fragment_counts <- SummarizedExperiment(assays = list(counts = my_counts_matrix),
                                          rowRanges = peaks)
  
  condition2 = all(colnames(my_counts_matrix) == sampleinfo[,SampleID])
  stopifnot(condition == "TRUE")
  
  sample_depth <- dplyr::select(sampleinfo, SampleDepth)
  fragment_counts[["depth"]] <- sample_depth[,1]
  fragment_counts[["depth"]]
  
  message("Adding GC Bias. it is required an input of a genome sequence to get background peaks based on the GC content and number of fragments across all samples  \n")
  fragment_counts <- addGCBias(fragment_counts, 
                               genome = BSgenome.Hsapiens.UCSC.hg38)
  
  
  message(paste("Filtering the samples and counts at minimum depth ",min_depth," minimum peaks ",min_peaks))
  counts_filtered <- filterSamples(fragment_counts,
                                   min_depth = min_depth,
                                   min_in_peaks = min_peaks,
                                   shiny = FALSE)
  
  filtering_plot <- filterSamplesPlot(fragment_counts,
                                      min_depth = min_depth,
                                      min_in_peaks = min_peaks, 
                                      use_plotly = FALSE)
  
  message("Saving the filtered plot\n")
  pdf(paste(saveDir,"filtering_plot.pdf",sep = ""),width = 5, height = 5)
  print(filtering_plot)
  dev.off()
  
  message("filtering out the samples\n")
  ix <- filterSamples(fragment_counts, 
                      min_depth = min_depth,
                      min_in_peaks = min_peaks,
                      ix_return = TRUE, shiny = FALSE)
  
  # Since counts_filtered has rownames so below line will not work. 
  # Check this if everything is correct we can move forward 
  # table(as.data.frame(getFragmentsPerPeak(counts_filtered)==0)) if all is false then move forward
  # without running the below line otherwise modify it
  # counts_filtered_remove_0 <- counts_filtered[-which(getFragmentsPerPeak(counts_filtered) == 0),] 
  
  message("fetching the details of the motifs for the filtered counts using Jaspar Motif\n")
  motifs <- getJasparMotifs()
  
  # The function matchMotifs from the motifmatchr package finds which peaks contain which motifs. By default, 
  # it returns a SummarizedExperiment object, which contains a sparse matrix (motif_ix@assays@data$motifMatches) indicating motif match or not.
  # The function requires an input of a genome sequence, which can be provided as a BSgenome, FaFile, or DNAStringSet object.
  # Can put the p.cutoff also for stringency,The default value is 0.00005, which tends to give reasonable numbers of motif matches for human motifs.
  
  # Instead of returning just motif matches, the function can also return additional matrices (stored as assays) with the number of motif matches
  # per peak and the maximum motif score per peak. For this additional information, use out = scores. To return the actual positions of motif
  # matches, use out = positions. Either the output with out = matches or out = scores can be passed to the computeDeviations function.
  message("Matching Motifs with the peaks\n")
  motif_ix <- matchMotifs(motifs, counts_filtered,
                          genome = BSgenome.Hsapiens.UCSC.hg38)
  
  message("saving motif RDS match annotation\n")
  saveRDS(motif_ix, paste(saveDir,"motif_ix.RDS",sep = ""))
  
  # The w paramter controls how similar background peaks should be. The bs parameter controls the precision with which the similarity is 
  # computed; increasing bs will make the function run slower. Sensible default parameters are chosen for both. It has 50 iterations
  message("Preparing the background peaks\n")
  bg <- getBackgroundPeaks(object = counts_filtered)
  
  
  # The function computeDeviations returns a SummarizedExperiment with two “assays”. The first matrix (accessible via deviations(dev) or 
  # assays(dev)$deviations) will give the bias corrected “deviation” in accessibility for each set of peaks (rows) for each cell or sample 
  # (columns). This metric represent how accessible the set of peaks is relative to the expectation based on equal chromatin accessibility 
  # profiles across cells/samples, normalized by a set of background peak sets matched for GC and average accessability. The second matrix 
  # (deviationScores(dev) or assays(deviations)$z) gives the deviation Z-score, which takes into account how likely such a score would occur if 
  # randomly sampling sets of peaks with similar GC content and average accessibility.
  
  message("Computing the devaition over the filtered peaks\n")
  dev_motif <- computeDeviations(object = counts_filtered, 
                                 annotations = motif_ix,
                                 background_peaks = bg)
  
  message("saving deviation motif RDS match annotation\n")
  saveRDS(dev_motif, paste(saveDir,"dev_motif.RDS",sep = ""))
  
  # The function computeVariability returns a data.frame that contains the variability (standard deviation of the z scores computed above across 
  # all cell/samples for a set of peaks), bootstrap confidence intervals for that variability (by resampling cells/samples), 
  # and a p-value for the variability being greater than the null hypothesis of 1.
  message("Computing the variability from the computed deviation\n")
  variability_motif <- computeVariability(dev_motif)
  
  message("saving variability motif RDS \n")
  saveRDS(variability_motif, paste(saveDir,"variability_motif.RDS",sep = ""))
  
  write.table(variability_motif, paste(saveDir,"variability_motif.txt",sep = ""), col.names = T, row.names = T, quote = F, sep = "\t")
  
  message("Plotting the variability plot\n")
  pdf(paste(saveDir,"variability_motif.pdf",sep = ""), width = 7, height = 7)
  print(plotVariability(variability_motif, use_plotly = FALSE))
  dev.off()
  
  stopifnot(all(sampleinfo$SampleID2 == colnames(dev_motif@assays@data$deviations))=="TRUE")
  dev_motif[["Group"]] <- sampleinfo[,Group]
  
  message("Making the tsne plots at 3 different threshold 1,1.5, and 2 \n")
  set.seed(2021)
  threshold = c(1,1.5,2)
  for (i in 1:length(threshold)) {
    tsne_results <- deviationsTsne(dev_motif, threshold = threshold[i], perplexity = 10)
    colnames(tsne_results) <- c("Dim1","Dim2")
    tsne_result_df <- as.data.frame(tsne_results)
    tsne_result_df$Group <- sampleinfo[,Group]
    tsne_result_df$Sample <- sampleinfo[,SampleID_2]
    p <- ggplot(tsne_result_df, 
                aes(Dim1,Dim2, label=Sample, color=Group)) +
      geom_point(size=2) + 
      geom_text_repel(size=2) + 
      ggtitle(paste("ChromVar Motif tsne threshold ",threshold[i], sep = "")) +
      theme(plot.title = element_text(hjust = 0.5), 
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(saveDir,"chromVar_tsne_group_Motif_threshold_",threshold[i],".pdf",sep = ""),width = 8, height=6)
    print(p)
    dev.off()
  }
  
  
  ## Doing the same thing with vierstra
  message("After Running the Jaspar Motif, we would do the same with the Vierstra Motifs")
  motifPWMs <- readRDS("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/Vierstra-Human-Motifs.rds") ## if not available download the Vierstra Motif RDS file
  
  message("Matching Vierstra Motifs with the peaks")
  motif_ix_2 <- matchMotifs(motifPWMs, counts_filtered, 
                            genome = BSgenome.Hsapiens.UCSC.hg38)
  
  message("saving motif RDS match annotation for Vierstra\n")
  saveRDS(motif_ix_2, paste(saveDir,"motif_ix_Vierstra.RDS",sep = ""))
  
  # if instead of knowing all the known motif and interested in motif of lenght 6 run below code
  # kmer_ix <- matchKmers(6, counts_filtered, 
  #                       genome = BSgenome.Hsapiens.UCSC.hg38)
  
  # We will get the deviation score and z score with compute deviations read the ChromVar green leaf if you want to know more
  message("Getting the background peaks")
  bg <- getBackgroundPeaks(object = counts_filtered) # background is same for both Jaspar and Vierstra
  
  message("Computing deviations")
  dev_vierstra <- computeDeviations(object = counts_filtered,
                                    annotations = motif_ix_2,
                                    background_peaks = bg)
  
  message("saving deviation motif RDS match annotation\n")
  saveRDS(dev_vierstra, paste(saveDir,"dev_vierstra.RDS",sep = ""))
  
  message("Computing variability")
  variability_vierstra <- computeVariability(dev_vierstra)
  
  message("saving variability motif RDS \n")
  saveRDS(variability_vierstra, paste(saveDir,"variability_vierstra.RDS",sep = ""))
  
  write.table(variability_vierstra, paste(saveDir,"variability_vierstra.txt",sep = ""), col.names = T, row.names = T, quote = F, sep = "\t")
  
  # deviations(dev_vierstra)
  # deviationScores(dev_vierstra)
  
  message("saving the Vierstra variability plot")
  pdf(paste(saveDir,"variability_Vierstra.pdf",sep = ""), width = 7, height = 7)
  print(plotVariability(variability_vierstra, use_plotly = FALSE))
  dev.off()
  
  stopifnot(all(sampleinfo$SampleID2 == colnames(dev_vierstra@assays@data$deviations))=="TRUE")
  dev_vierstra[["Group"]] <- sampleinfo[,Group]
  
  message("Making the Vierstra tsne plots at 3 different threshold 1,1.5, and 2 \n")
  set.seed(2021)
  threshold = c(1,1.5,2)
  for (i in 1:length(threshold)) {
    tsne_results <- deviationsTsne(dev_vierstra, threshold = threshold[i], perplexity = 10)
    colnames(tsne_results) <- c("Dim1","Dim2")
    tsne_result_df <- as.data.frame(tsne_results)
    tsne_result_df$Group <- sampleinfo[,Group]
    tsne_result_df$Sample <- sampleinfo[,SampleID_2]
    p <- ggplot(tsne_result_df, 
                aes(Dim1,Dim2, label=Sample, color=Group)) +
      geom_point(size=2) + 
      geom_text_repel(size=2) + 
      ggtitle(paste("ChromVar Vierstra tsne threshold ",threshold[i], sep = "")) +
      theme(plot.title = element_text(hjust = 0.5), 
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey"))
    
    pdf(paste(saveDir,"chromVar_tsne_group_Vierstra_threshold_",threshold[i],".pdf",sep = ""),width = 8, height=6)
    print(p)
    dev.off()
  }
  
  
  
  ## Making the heatmap
  top_50_variability <- variability_motif[order(-variability_motif$variability),] %>% head(50)
  deviation_score <- dev_motif@assays@data$z
  
  Group2 <- unique(sampleinfo[,Group])
  for (i in 1:length(Group2)){
    name <- paste(Group2[i],"index", sep = "_")
    index <- grep(Group2[i], sampleinfo[,Group])
    assign(name, index)
  }
  
  Group2_index <- paste(Group2,"index", sep = "_")
  
  # making a median rlog for each group according to the target Group2 
  message("Making an empty dataframe\n")
  deviation_score.meds <- data.frame(matrix(ncol = length(Group2_index), nrow = nrow(deviation_score))) # making an empty dataframe
  rownames(deviation_score.meds) <- rownames(deviation_score)
  colnames(deviation_score.meds) <- Group2
  i=0
  j=0
  
  message("Calculating the median for the group\n")
  for (i in 1:nrow(deviation_score.meds)){
    for(j in 1:length(Group2_index)){
      deviation_score.meds[i,j] <- median(as.matrix(deviation_score[i,get(Group2_index[j])]))
    }
  }
  
  message("Saving the median matrix\n")
  write.table(deviation_score.meds,
              file = paste(saveDir,"deviation_score_median_jaspar_motif.txt",sep = ""),
              sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  
  
  write.table(dev_motif@assays@data$deviations,
              file = paste(saveDir,"deviation_score_jaspar_motif.txt",sep = ""),
              sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  
  
  message("Matching the top most variable motifs in the Jaspar\n")
  deviation_score_top_50 <- deviation_score.meds[match(rownames(top_50_variability),rownames(deviation_score.meds)),]
  
  deviation_score_top_50_scaled <- scale(deviation_score_top_50)
  
  # column_anno <- HeatmapAnnotation(Group = targets_old_and_young_vzv$Group,
  #                                  Age = targets_old_and_young$Age,
  #                                  col = list(Group = c("old" = "darksalmon", 
  #                                                       "young" = "lightblue")))
  h <- Heatmap(deviation_score_top_50_scaled, 
               name = "median dev",
               column_title = "samples",
               row_title = "TFs",
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               cluster_columns = TRUE,
               cluster_rows = TRUE
               #top_annotation = column_anno
  )
  
  pdf(paste(saveDir,"deviation_Score_median_heatmap_motif_top50_Jaspar_clustered.pdf",sep = ""), height = 7, width = 6)
  print(h)
  dev.off()
  
  # ordering the deviation for top 50 scaled naive, cm, em, ebv
  deviation_score_top_50_scaled_order <- deviation_score_top_50_scaled[,col_arrange]
  h <- Heatmap(deviation_score_top_50_scaled_order, 
               name = "median dev",
               column_title = "samples",
               row_title = "TFs",
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               cluster_columns = FALSE,
               cluster_rows = TRUE
               #top_annotation = column_anno
  )
  
  pdf(paste(saveDir,"deviation_Score_median_heatmap_motif_top50_Jaspar_unclustered.pdf",sep = ""), height = 7, width = 6)
  print(h)
  dev.off()
  
  ### Now vierstra
  top_50_variability <- variability_vierstra[order(-variability_vierstra$variability),] %>% head(50)
  deviation_score <- dev_vierstra@assays@data$deviations
  
  Group2 <- unique(sampleinfo[,Group])
  for (i in 1:length(Group2)){
    name <- paste(Group2[i],"index", sep = "_")
    index <- grep(Group2[i], sampleinfo[,Group])
    assign(name, index)
  }
  
  Group2_index <- paste(Group2,"index", sep = "_")
  
  # making a median rlog for each group according to the target Group2 
  deviation_score_vierstra.meds <- data.frame(matrix(ncol = length(Group2_index), nrow = nrow(deviation_score))) # making an empty dataframe
  rownames(deviation_score_vierstra.meds) <- rownames(deviation_score)
  colnames(deviation_score_vierstra.meds) <- Group2
  i=0
  j=0
  
  for (i in 1:nrow(deviation_score_vierstra.meds)){
    for(j in 1:length(Group2_index)){
      deviation_score_vierstra.meds[i,j] <- median(as.matrix(deviation_score[i,get(Group2_index[j])]))
    }
  }
  
  write.table(deviation_score_vierstra.meds, 
              file = paste(saveDir,"deviation_score_median_vierstra_motif.txt",sep = ""),
              sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  
  write.table(dev_vierstra@assays@data$deviations,
              file = paste(saveDir,"deviation_score_vierstra_motif.txt",sep = ""),
              sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  
  deviation_score_top_50 <- deviation_score_vierstra.meds[match(rownames(top_50_variability),rownames(deviation_score_vierstra.meds)),]
  deviation_score_top_50_scaled <- scale(deviation_score_top_50)
  
  h <- Heatmap(deviation_score_top_50_scaled,
               name = "dev scaled",
               column_title = "samples",
               row_title = "TFs",
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               cluster_columns = TRUE,
               cluster_rows = TRUE
  )
  
  pdf(paste(saveDir,"deviation_Score_heatmap_vierstra_top50_scaled_clustered.pdf",sep = ""), height = 8, width = 7)
  print(h)
  dev.off()
  
  deviation_score_top_50_scaled_order <- deviation_score_top_50_scaled[,col_arrange]
  h <- Heatmap(deviation_score_top_50_scaled_order, 
               name = "median dev",
               column_title = "samples",
               row_title = "TFs",
               row_names_gp = gpar(fontsize = 8),
               column_names_gp = gpar(fontsize = 8),
               cluster_columns = FALSE,
               cluster_rows = TRUE
               #top_annotation = column_anno
  )
  
  pdf(paste(saveDir,"deviation_Score_heatmap_vierstra_top50_scaled_unclustered.pdf",sep = ""), height = 8, width = 7)
  print(h)
  dev.off()
  
  message("Saving the ChromVar session")
  save.session(paste(saveDir,"ChromVar_vzv.Rda",sep = ""))
}







