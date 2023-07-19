##### Running Homer ########
### Running for comparison 1 
# final_2 is the kmeans generated during the kmeans(df_scale, clus_num, nstart = 25) 
# df_scale is the cnts2.norm.sig.z.meds
# colindex_arrange could be the column the way you want to arrange in the heatmap for the df_scale
homer_ATAC <- function(saveDir,final_2,df_scale,colindex_arrange,column_anno=NULL){
  ## Making a heatmap for the optimal cluster of peaks in rows and groups in columns
  peaks_arranged <- as.data.frame(final_2$cluster[order(as.data.frame(final_2$cluster))])
  colnames(peaks_arranged) <- "clusters"
  peaks_arranged$clusters <- paste("C",peaks_arranged$clusters,sep="")
  table(peaks_arranged$clusters)
  
  df_scale_ordered <- df_scale[match(rownames(peaks_arranged), rownames(df_scale)),]
  
  ## The marking and labels for the row anno you will get after you run the homer below there is the code
  ## so first run the homer then you will get both these values then we can run this
  cluster_and_freq <- as.data.frame(table(peaks_arranged$clusters))
  sorted_clusters <- paste("C",sort(as.numeric(gsub("C","",cluster_and_freq$Var1))), sep = "")
  match(sorted_clusters, cluster_and_freq$Var1)
  cluster_and_freq <- cluster_and_freq[match(sorted_clusters, cluster_and_freq$Var1),]
  frequency <- cluster_and_freq$Freq
  rm(marking)
  marking <- vector()
  marking[1] <- (frequency[1]/2)
  for (i in 2:length(frequency)){print(i)
    marking[i] <- (sum(frequency[1:i-1]) + (frequency[i]/2))
  }
  marking <- round(marking)
  
  # label you can get from the below code that makes multiple plot for Homer it is running with the for loop
  
  ### Prioritizing the TFs according to the cluster from the Homer result.
  files <- vector()
  for(i in 1:length(unique(peaks_arranged$clusters))){
    files[i] = list.files(path = paste(saveDir,"Homer/","C",i,"_coordinates_df",sep = ""),
                          pattern = "knownResults.txt",
                          recursive = T,
                          full.names = T)
  }
  
  setAs("character", "myDouble", function(from)as.double(sub('^(-?) e','\\13e',from))) # to read the e values in write.delim
  
  i=0
  rm(plot_list_4, plot_list_5, plot_list_6, label_vec)
  plot_list_4 <- list()
  plot_list_5 <- list()
  plot_list_6 <- list()
  label_vec <- vector() # making this vector for above heatmap that include the TFs.
  
  for (i in 1:length(files)){
    homer_tb <- read.delim(files[i], sep = "\t", header = T,
                           colClasses=list ("character", "character", "myDouble","myDouble",
                                            "myDouble", "double", "character","double","character"))
    colnames(homer_tb) <- c("Motif", "Consensus", "Pvalue", "logPvalue", "qvalue", "foreground",
                            "foreground_percent", "background", "background_percent")
    homer_tb$Enrichment <- homer_tb$foreground/homer_tb$background
    homer_tb$Motif <- gsub("?","QuestionMark",homer_tb$Motif, fixed = TRUE)
    rm(TF)
    TF <- vector()
    j=0
    for(j in 1:nrow(homer_tb)){
      TF[j] <- gsub("[\\(\\)]", "", regmatches(homer_tb$Motif, gregexpr("\\(.*?\\)", homer_tb$Motif))[[j]])
    }
    TF_uniq <- unique(TF)
    rm(idx)
    idx <- vector()
    for (k in 1:length(TF_uniq)) {
      TF_uniq_exact <- paste("^",TF_uniq[k],"$", sep = "")
      idx[k] <- grep(TF_uniq_exact,TF)[1]
    }
    homer_tb_2 <- homer_tb[idx,]
    homer_tb_2$TF_family <- TF_uniq
    homer_tb_2 <- remove.factors(homer_tb_2)
    homer_tb_2$TF_and_family <- gsub("\\/.*","",homer_tb_2$Motif)
    #homer_tb_2$TF_family <- factor(homer_tb_2$TF_family, levels = homer_tb_2$TF_family[order(homer_tb_2$foreground)])
    homer_tb_2$TF_and_family <- factor(homer_tb_2$TF_and_family, levels = homer_tb_2$TF_and_family[order(homer_tb_2$foreground)])
    p <- ggplot(homer_tb_2, aes(y=TF_and_family, x=foreground, color = logPvalue, size=Enrichment)) + 
      geom_point() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y=element_text(size=9),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey")) +
      ggtitle(paste("C",i," Homer Result",sep = ""))
    plot_list_4[[i]] <- p
    
    homer_tb_2 <- remove.factors(homer_tb_2)
    homer_tb_2$TF_and_family <- factor(homer_tb_2$TF_and_family, levels = homer_tb_2$TF_and_family[order(homer_tb_2$Enrichment)])
    #homer_tb_2$TF_family <- factor(homer_tb_2$TF_family, levels = homer_tb_2$TF_family[order(homer_tb_2$Enrichment)])
    q <- ggplot(homer_tb_2, aes(y=TF_and_family, x=Enrichment, color = logPvalue, size=foreground)) + 
      geom_point() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y=element_text(size=9),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey")) +
      ggtitle(paste("C",i," Homer Result",sep = ""))
    plot_list_5[[i]] <- q
    
    homer_tb_2 <- remove.factors(homer_tb_2)
    #homer_tb_2$TF_family <- factor(homer_tb_2$TF_family, levels = homer_tb_2$TF_family[order(-homer_tb_2$logPvalue)])
    homer_tb_2$TF_and_family <- factor(homer_tb_2$TF_and_family, levels = homer_tb_2$TF_and_family[order(homer_tb_2$logPvalue)])
    l <- ggplot(homer_tb_2, aes(y=TF_and_family, x=-logPvalue, color = Enrichment, size=foreground)) + 
      geom_point() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y=element_text(size=9),
            panel.background = element_rect(fill = 'white', colour = 'black'),
            panel.grid.minor = element_line(colour = "grey"),
            panel.grid.major = element_line(colour = "grey")) +
      ggtitle(paste("C",i," Homer Result",sep = ""))
    plot_list_6[[i]] <- l
    
    # remove the factor for the pvalues as well now we need the top 3 values
    homer_tb_2 <- remove.factors(homer_tb_2)
    #homer_tb_2$TF_family <- factor(homer_tb_2$TF_family, levels = homer_tb_2$TF_family[order(homer_tb_2$logPvalue)])
    homer_tb_2$TF_and_family <- factor(homer_tb_2$TF_and_family, levels = homer_tb_2$TF_and_family[order(homer_tb_2$logPvalue)])
    selected <- homer_tb_2[homer_tb_2$Pvalue < 0.05,"TF_and_family"]
    if(length(selected)>3){
      label_vec[i] <- paste(selected[1:3], collapse = ";")
    }
    else{
      label_vec[i] <- paste(selected,collapse = ";")
    }
    #label_vec[i] <- paste(levels(homer_tb_2$TF_family)[1:3], collapse = ";")
    rm(TF_uniq_exact, TF_uniq, TF)
  }

  pdf(paste(saveDir,"Homer/","Homer_TF_Foreground_order.pdf",sep = ""), width = 20, height = 30)
  require(gridExtra)
  grid.arrange(plot_list_4[[1]],plot_list_4[[2]],plot_list_4[[3]],plot_list_4[[4]],plot_list_4[[5]],
               plot_list_4[[6]],plot_list_4[[7]],
               plot_list_4[[8]],plot_list_4[[9]],plot_list_4[[10]],plot_list_4[[11]],
               plot_list_4[[12]],plot_list_4[[13]],plot_list_4[[14]],
               plot_list_4[[15]],plot_list_4[[16]],plot_list_4[[17]],
               plot_list_4[[18]],plot_list_4[[19]],plot_list_4[[20]],plot_list_4[[21]],plot_list_4[[22]],plot_list_4[[23]],
               plot_list_4[[24]],plot_list_4[[25]],plot_list_4[[26]],plot_list_4[[27]],plot_list_4[[28]],plot_list_4[[29]],
               ncol=3)
  dev.off()
  
  pdf(paste(saveDir,"Homer/","Homer_TF_Enrichment_order.pdf",sep = ""), width = 20, height = 30)
  require(gridExtra)
  grid.arrange(plot_list_5[[1]],plot_list_5[[2]],plot_list_5[[3]],plot_list_5[[4]],plot_list_5[[5]],
               plot_list_5[[6]],plot_list_5[[7]],
               plot_list_5[[8]],plot_list_5[[9]],plot_list_5[[10]],plot_list_5[[11]],
               plot_list_5[[12]],plot_list_5[[13]],plot_list_5[[14]],
               plot_list_5[[15]],plot_list_5[[16]],plot_list_5[[17]],
               plot_list_5[[18]],plot_list_5[[19]],plot_list_5[[20]],plot_list_5[[21]],plot_list_5[[22]],plot_list_5[[23]],
               plot_list_5[[24]],plot_list_5[[25]],plot_list_5[[26]],plot_list_5[[27]],plot_list_5[[28]],plot_list_5[[29]],
               ncol=3)
  dev.off()
  
  
  pdf(paste(saveDir,"Homer/","Homer_TF_pvalue_order.pdf",sep = ""), width = 20, height = 30)
  require(gridExtra)
  grid.arrange(plot_list_6[[1]],plot_list_6[[2]],plot_list_6[[3]],plot_list_6[[4]],plot_list_6[[5]],
               plot_list_6[[6]],plot_list_6[[7]],
               plot_list_6[[8]],plot_list_6[[9]],plot_list_6[[10]],plot_list_6[[11]],
               plot_list_6[[12]],plot_list_6[[13]],plot_list_6[[14]],
               plot_list_6[[15]],plot_list_6[[16]],plot_list_6[[17]],
               plot_list_6[[18]],plot_list_6[[19]],plot_list_6[[20]],plot_list_6[[21]],plot_list_6[[22]],plot_list_6[[23]],
               plot_list_6[[24]],plot_list_6[[25]],plot_list_6[[26]],plot_list_6[[27]],plot_list_6[[28]],plot_list_6[[29]],
               ncol=3)
  dev.off()
  
  
  
  row_anno <- rowAnnotation(Clusters =  peaks_arranged$clusters,
                            col = list(Clusters = c("C1" = "aquamarine3", 
                                                    "C2" = "brown2",
                                                    "C3" = "chocolate1",
                                                    "C4" = "dodgerblue1",
                                                    "C5" = "lightcoral",
                                                    "C6" = "lightgreen",
                                                    "C7" = "burlywood4",
                                                    "C8" = "firebrick3",
                                                    "C9" = "cyan4",
                                                    "C10" = "chartreuse3",
                                                    "C11" = "bisque3",
                                                    "C12" = "darksalmon",
                                                    "C13" = "darkseagreen4",
                                                    "C14" = "mediumorchid3",
                                                    "C15" = "lightslategray",
                                                    "C16" = "aquamarine2", 
                                                    "C17" = "brown1",
                                                    "C18" = "chocolate2",
                                                    "C19" = "dodgerblue2",
                                                    "C20" = "coral",
                                                    "C21" = "green",
                                                    "C22" = "burlywood3",
                                                    "C23" = "firebrick2",
                                                    "C24" = "cyan3",
                                                    "C25" = "chartreuse2",
                                                    "C26" = "bisque2",
                                                    "C27" = "salmon",
                                                    "C28" = "darkseagreen2",
                                                    "C29" = "mediumorchid2"
                            )),
                            month = anno_mark(at= marking,labels = label_vec))
  
  h=Heatmap(df_scale_ordered,
            name = "scaled", #title of legend
            column_title = "Groups", row_title = "Peaks",
            row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
            column_names_gp = gpar(fontsize = 10),
            cluster_columns = FALSE,
            cluster_rows = FALSE,
            #rect_gp = gpar(col = "white", lwd = 1),
            right_annotation = row_anno,
            top_annotation = column_anno
            #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
  )
  
  pdf(paste(saveDir,"Homer/","Homer_cluster_top_TFs_clustered.pdf",sep = ""), width = 15, height = 10)
  print(h)
  dev.off()
  
  p=Heatmap(df_scale_ordered,
            name = "scaled", #title of legend
            column_title = "Groups", row_title = "Peaks",
            row_names_gp = gpar(fontsize = 0.0001), # Text size for row names
            column_names_gp = gpar(fontsize = 10),
            cluster_columns = TRUE,
            cluster_rows = FALSE,
            #rect_gp = gpar(col = "white", lwd = 1),
            right_annotation = row_anno,
            top_annotation = column_anno,
            column_split = factor(c(rep("A", 19), rep("B",19), rep("C",18)), levels = c("A","B","C")),
            #column_split = factor(c(rep("A", 9), rep("B",8)), levels = c("A","B")),
  )
  
  
  
  
  pdf(paste(saveDir,"Homer/","Homer_cluster_top_TFs_splitted.pdf",sep = ""), width = 15, height = 10)
  print(p)
  dev.off()
  
}

  
