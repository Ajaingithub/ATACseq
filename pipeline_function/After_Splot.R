# This function is run after running the splot_ATACseq.R i.e. after we have identified the number of overlaps that is required for filtering the peaks
# This function is used for filtering the peaks so that the matrix could be further use for the downstream analysis
# Args:
# All_merge_path: The path of the all_merge.txt (output for splot_ATACseq.R)
# overlaps: Depending on the Splot generated from splot_ATAC.R decide the number of overlaps to remove the reads with less than overlap numnber
# saveDir: save directory path
# bam_path: sorted remove duplicate bam files for all the samples present in the all_merge.txt
# blacklist_region: To remove the blacklisted region from the ATAC data
# sampleinfo: samples metadata 

After_Splot <- function(All_merge_path,overlaps,saveDir,bam_path,blacklist_region,sampleinfo){
  message("loading All merge file")
  All_merge <- read.table(All_merge_path, sep = "\t", header = TRUE)
  
  message("Filtering the reads according to the splot\n")
  All_merge_filtered<-All_merge[!All_merge$n.overlaps %in% 1:overlaps,] # This removes all those region with less than 6 overlaps
  All_combined_merged<-bedr.merge.region(All_merge_filtered,check.chr = FALSE) # Takes time
  
  # Taking only autosomal chromosome since mitochondrial DNA is opened, more Tn5 transposons are
  # attached to it as it is more accessible. It may contain 34% of the mitochondrial DNA. there 
  # are no ATAC-seq peaks of interest in the mitochondrial genome, these reads will only complicate the subsequent steps.
  #-no-discordant --no-dovetail --no-mixed --no-unal
  # Also sex chromosome will also bias the analysis so they are also been removed from the sample
  message("Removing the sex chromosome and chrM\n")
  All_combined_merged <- All_combined_merged[All_combined_merged$chr %in% c("chr1", "chr2", "chr3","chr4", "chr5", "chr6",
                                                                            "chr7","chr8","chr9","chr10","chr11","chr12",
                                                                            "chr13","chr14","chr15","chr16","chr17","chr18",
                                                                            "chr19","chr20","chr21","chr22"),]
  
  #### Read counts mapping ####
  # con <- file(paste(saveDir,"Rsubread_",RN,"_Subset.log",sep = ""))
  # sink(con,append = TRUE,split = TRUE)
  
  # Rsubread package Alignment, quantification and analysis of RNA sequencing data (including both bulk RNA-seq and scRNA-seq)
  # and DNA sequenicng data (including ATAC-seq, ChIP-seq, WGS, WES etc). Includes functionality for read
  # mapping, read counting, SNP calling, structural variant detection and gene fusion discovery. Can be 
  # applied to all major sequencing techologies and to both short and long sequence reads.
  
  # rtracklayer Extensible framework for interacting with multiple genome browsers (currently UCSC built-in)
  #and manipulating annotation tracks in various formats (currently GFF, BED, bedGraph, BED15, WIG, BigWig 
  #and 2bit built-in). The user may export/import tracks to/from the supported browsers, as well as query 
  #and modify the browser state, such as the current viewport.
  
  # When working on large lists or data.frames, it might be both time and memory consuming to convert them to a data.table using
  # as.data.table(.), as this will make a complete copy of the input object before to convert it to a data.table. The setDT function
  # takes care of this issue by allowing to convert lists - both named and unnamed lists and data.frames by reference instead. That
  # is, the input object is modified in place, no copy is being made. 
  setDT(All_combined_merged,keep.rownames=TRUE)[]
  colnames(All_combined_merged)[1]<-"GeneID" 
  All_combined_merged$strand <- rep("-",nrow(All_combined_merged))
  All_combined_merged_for_featureCounts<-All_combined_merged[,c(1,2,3,4,6)] #removed the names column
  
  #It provides the number of reads mapping to each genomic feature, 
  #Note that annot.ext will override annot.inbuilt if both provided.
  message("Performing Feature Counts \n")
  All_merge_counts <-featureCounts(files = bam_path, annot.ext = All_combined_merged_for_featureCounts, 
                                   useMetaFeatures = FALSE, isPairedEnd = TRUE, nthreads = 50) 
  #The argument useMetaFeatures specifies whether read summarization should be performed at feature level 
  #or at meta-feature level. A feature represents a continuous genomic region and a meta-feature is a group 
  #of features. For instance, an exon is a feature and a gene comprising one or more exons is a meta-feature. 
  
  write.table(All_merge_counts$counts, file = paste(saveDir,"Called_peaks_table_all.txt",sep = "/"),
              sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
  write.table(All_merge_counts$stat, file = paste(saveDir,"Called_peaks_table_all_stat.txt",sep = "/"),
              sep = "\t",col.names = TRUE,row.names = TRUE,quote = FALSE)
  
  ### Create Tracks ###
  write.table(t(All_merge_counts$stat[1,2:ncol(All_merge_counts$stat)]),file = paste(saveDir,"SampleDepth.txt", sep = "/"),sep = "\t",col.names = FALSE,row.names = FALSE,quote = FALSE)

  cnts_orig <- read.table(paste(saveDir,"Called_peaks_table_all.txt",sep = ""))
  cnts <- read.table(paste(saveDir,"Called_peaks_table_all.txt", sep = ""))
  colnames(cnts)<-sampleinfo$SampleID2       ### Reconfirm the nomenclature!!!
  targets_sorted<-sampleinfo
  
  library(data.table)
  cnts_sorted<-cnts
  setDT(cnts_sorted)
  setcolorder(cnts_sorted, as.character(targets_sorted$SampleID2))
  
  cnts_sorted<-as.data.frame(cnts_sorted)
  print(head(targets_sorted$SampleID2))
  print(head(colnames(cnts_sorted)))
  
  ### Remove BlackList Regions #####
  cnts <- read.table(paste(saveDir,"Called_peaks_table_all.txt", sep = ""))
  row.names(cnts_sorted)<-row.names(cnts)
  cnts.bed<-row.names(cnts)
  cnts.bed<-as.data.frame(cnts.bed)
  cnts.bed[,c("V1","V2")]<-str_split_fixed(cnts.bed$cnts.bed, ":", 2)
  cnts.bed[,c("V2","V3")]<-str_split_fixed(cnts.bed$V2, "-", 2)
  cnts.bed<-cnts.bed[,c("V1","V2","V3")]
  write.table(x = cnts.bed,file = paste(saveDir,"cnts.bed", sep = ""),sep = "\t",row.names = F,col.names = F,quote = F)
  
  # If this file is absent /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run6_and_Run7/wgEncodeHg38ConsensusSignalArtifactRegions.bed then run this 
  # cd /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Rohit/VirusSpecificTcells/CD4/Huimin_Rerun_Sept21/Run6/Analysis/Downstream
  # wget https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
  # gunzip ENCFF356LFX.bed.gz
  # mv ENCFF356LFX.bed wgEncodeHg38ConsensusSignalArtifactRegions.bed
  # module load bedtools/2.27.1
  # permanently added in the mayobiotools
  
  setwd(saveDir)
  system(paste("bedtools subtract -a ./cnts.bed -b ",blacklist_region," > ./cnts_noblacklist.bed",sep = ""))
  
  cnts_noblacklist<-read.table(paste(saveDir,"cnts_noblacklist.bed",sep = ""),header = F)
  cnts_noblacklist$names<-paste(cnts_noblacklist$V1, cnts_noblacklist$V2, sep=":")
  cnts_noblacklist$names<-paste(cnts_noblacklist$names, cnts_noblacklist$V3, sep="-")
  cnts_noblacklist_setdiff<-setdiff(rownames(cnts_sorted),cnts_noblacklist$names) # different rownames
  cnts_noblacklist_data<-cnts_sorted[rownames(cnts_sorted) %in% cnts_noblacklist$names, ] # removing the blacklisted region
  write.table(x = cnts_noblacklist_data,file = paste(saveDir,"cnts_noblacklist_data.txt",sep = ""),
              sep = "\t",row.names = T,col.names = T,quote = F)
  write.table(x = as.matrix(colSums(cnts_noblacklist_data)), file = paste(saveDir,"Sample_Names_Depth_NoBL.txt", sep = ""),
              sep = "\t",row.names = T,col.names = F,quote = F)
}

