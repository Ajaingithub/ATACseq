splot <- function(peak_path, build = "hg38", species = "human", saveDir){
  message("Reading all the narrow peaks\n")
  All <- lapply(peak_path,function(i){
    read.delim(i, header = FALSE, sep = "\t")
  })
  
  peak_path <- basename(peak_path)
  peak_path <- gsub("-",".",peak_path)
  names(All) <- gsub("_peaks.narrowPeak","",peak_path)
  
  message("Making a bed file from the narrow peaks\n")
  for (i in 1:length(All)){
    All[[i]][1]<-lapply(All[[i]][1],as.character)
    All[[i]]<-All[[i]][,1:3]
    All[[i]]<-convert2bed(All[[i]],check.chr = FALSE)
  }
  
  ## Joined all the bed files
  message("Joining all the bed files, This takes time\n")
  All_merge<-bedr.join.multiple.region(x = All,fraction.overlap = 0.5, 
                                       species = species, 
                                       build = build,
                                       check.chr = FALSE)
  rm(All)
  All_Peak_Stats <- data.frame(Overlaps=unique(All_merge$n.overlaps))
  All_Peak_Stats$Peaks<-rep("NA",nrow(All_Peak_Stats)) # Just added NA to the Peak column
  
  message("Adding the peaks for each overlap\n")
  for(i in 1:nrow(All_Peak_Stats)){
    temp<-All_merge[All_merge$n.overlaps %in% All_Peak_Stats$Overlaps[i],]
    All_Peak_Stats$Peaks[i]<-nrow(temp)
  }
  
  All_Peak_Stats$Overlaps <- as.numeric(as.character(All_Peak_Stats$Overlaps))
  All_Peak_Stats_sorted <- All_Peak_Stats[order(All_Peak_Stats$Overlaps),] 
  
  message("Saving the sorted peaks\n")
  write.table(All_Peak_Stats_sorted,file = paste(saveDir,"All_Peak_Stats_sorted.txt",sep = ""),sep = "\t",col.names = TRUE,row.names = FALSE,quote = FALSE)
  
  message("Saving the All merge peaks\n")
  write.table(All_merge,file = paste(saveDir,"All_merge.txt",sep = ""),sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
  
  All_Peak_Stats_sorted <- read.table(paste(saveDir,"All_Peak_Stats_sorted.txt", sep = ""), header = TRUE)
  
  message("Making an S plot to check for the overlaps\n")
  p <- ggplot(All_Peak_Stats_sorted, aes(Overlaps, Peaks)) + geom_line() + 
    theme_bw() + scale_x_continuous(breaks= seq(0,length(peak_path),by=2)) + 
    scale_y_continuous(breaks= seq(0,All_Peak_Stats_sorted[1,2]+10000,by=10000))
  
  pdf(paste(saveDir,"Peaks_Stats_S_plot.pdf",sep = ""), width = 8, height = 8)
  print(p)
  dev.off()
  
  message("Please check the splot at this location ",saveDir," and decide how many needs to be kept")
}

