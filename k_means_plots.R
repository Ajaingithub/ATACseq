## k_means heatmap
k_means_plot <- function(cluster_gap, df_scale, saveDir,clus_num,deseq_dataset){
  `%notin%` <- Negate(`%in%`)
  dir.create(paste(saveDir,"Clustering",sep = ""),showWarnings = FALSE)
  pdf(paste(saveDir,"Clustering/cluster_number.pdf",sep = ""), width = 7, height = 7)
  print(fviz_gap_stat(cluster_gap))
  dev.off()
  
  dds <- deseq_dataset
  
  set.seed(123)
  final_2 <- kmeans(df_scale, clus_num, nstart = 25) # based on the optimal clusters
  print(final_2)
  
  saveRDS(final_2,paste(saveDir,"Clustering/final_2.Rds",sep = ""))
  
  peaks_arranged <- as.data.frame(final_2$cluster[order(as.data.frame(final_2$cluster))])
  colnames(peaks_arranged) <- "clusters"
  peaks_arranged$clusters <- paste("C",peaks_arranged$clusters,sep="")
  peak_distribution <- table(peaks_arranged$clusters)
  
  print(peak_distribution)
  write.table(peak_distribution,
              paste(saveDir,"numbers/peaks_distribution.txt",sep = ""), sep = "\t",
              row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  # prepapring data for Homer. For each Homer run each cluster with sig region and all other consider as non sig region
  peaks_arranged$coordinated <- rownames(peaks_arranged)
  clusters <- unique(peaks_arranged$clusters)
  i=0
  for(i in 1:length(clusters)){
    name <- paste(clusters[i],"coordinates_df",sep="_")
    cluster_coordinate <- peaks_arranged[grep(paste(clusters[i],"$",sep = ""),peaks_arranged$clusters),"coordinated"]
    assign(name,cluster_coordinate)
    backgroundname <- paste(clusters[i],"_cluster_background",sep="")
    background_peaks <- rownames(assay(dds))[which(rownames(assay(dds)) %notin% cluster_coordinate)]
    assign(backgroundname,background_peaks)
  }
  
  coordinate_files <- ls(pattern="coordinates_df")
  for(i in 1:length(coordinate_files)){
    filename <- paste(saveDir,"Clustering/",coordinate_files[i],".txt",sep = "")
    write.table(get(coordinate_files[i]), filename,
                row.names = FALSE, col.names = FALSE,
                sep = "\t", quote = FALSE)
  }
  background_files <- ls(pattern = "_cluster_background")
  for (i in 1:length(background_files)) {
    filename <- paste(saveDir,"Clustering/",background_files[i],".txt",sep = "")
    write.table(get(background_files[i]), filename,
                row.names = FALSE, col.names = FALSE,
                sep = "\t", quote = FALSE)
  }
  message(paste("Please check this location ",saveDir,"Clustering/ for the plots and significant and background files for running Homer" ,sep = ""))
}









