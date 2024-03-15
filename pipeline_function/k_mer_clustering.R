### Runnung K-mer clustering for all #####
# This function will generate the gap statistics to identify how many cluster (K) can be generated
# Args: 
# Differential_peaks: Differential Peaks generated during the previous process differential_ATAC.R
# rlog_dds: rlog Normalized counts data
# sampleinfo: samples metadata
# Group: Column from the count matrix need to combine
# saveDir: save directory

kmeans_clustering_ATAC <- function(differential_peaks,rlog_dds,sampleinfo,Group,saveDir){
  cnts2.norm <- assay(rlog_dds)
  cnts2.norm.sig <- cnts2.norm[match(differential_peaks,rownames(cnts2.norm)),]
  cnts2.norm.sig.z <- t(apply(cnts2.norm.sig, 1, function(i) (i - mean(i))/sd(i)))
  write.table(cnts2.norm.sig.z,
              file = paste(saveDir,"Differential/significant_normalized_all_samples.txt",sep = ""),
              sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  #cnts2.norm.sig.z.meds <- t(apply(cnts2.norm.sig.z, 1, function(i) sapply(levels(targets$Group2), function(l) median(i[targets$Group2 == l]))))
  
  Group2 <- unique(sampleinfo[Group])[,]
  for (i in 1:length(Group2)){
    name <- paste(Group2[i],"index", sep = "_")
    index <- grep(Group2[i], sampleinfo[Group][,])
    assign(name, index)
  }
  
  Group2_index <- paste(Group2,"index", sep = "_")
  
  # making a median rlog for each group according to the target Group2 
  cnts2.norm.sig.z.meds <- data.frame(matrix(ncol = length(Group2_index), nrow = nrow(cnts2.norm.sig.z))) # making an empty dataframe
  rownames(cnts2.norm.sig.z.meds) <- rownames(cnts2.norm.sig)
  colnames(cnts2.norm.sig.z.meds) <- Group2
  i=0
  j=0
  
  for(i in 1:nrow(cnts2.norm.sig.z)){
    for(j in 1:length(Group2_index)){
      cnts2.norm.sig.z.meds[i,j] <- median(as.matrix(cnts2.norm.sig.z[i,get(Group2_index[j])]))
    }
  }
  
  #cnts2.norm.sig.z.meds <- t(apply(cnts2.norm.sig.meds, 1, function(i) (i - mean(i))/sd(i)))
  
  write.table(cnts2.norm.sig.z.meds,
              file = paste(saveDir,"Differential/median_counts_for_Group_significant_normalized.txt",sep = ""),
              sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  
  df_scale <- cnts2.norm.sig.z.meds
  #df_scale_ordered_2 <- df_scale[,c(3,1,2,4,5,8,6,7,9,10)] ##Arranged naive, cm, em, flu, ebv, vzv
  

  distance <- get_dist(t(df_scale))
  pdf(paste(saveDir,"Differential/distance_heatmap.pdf",sep = ""), width = 7, height = 7)
  print(fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07")))
  dev.off()
  
  set.seed(123)
  gap_stat_2 <- clusGap(df_scale, FUN = kmeans, nstart = 25, K.max = 20) #df_scale wonot get transposed as we are finding differential based on the peaks
  print(gap_stat_2)
  
  saveRDS(gap_stat_2,paste(saveDir,"Differential/cluster_gap.RDS",sep = ""))
}
