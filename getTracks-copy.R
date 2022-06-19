library(Gviz)
library(rtracklayer)
library(trackViewer)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Rohit/GCP/API/atacTracksV2.R")

Track_list<-list(
  Old_CM = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/memory_bw_transfer/old_cm_Average_memory.bw",
  Young_CM = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/memory_bw_transfer/young_cm_Average_memory.bw",
  bed = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/memory_bw_transfer/cm_old_vs_young_forground_down.bed"
)

Track_cols<-c("blue","green","yellow")

Genes_names<-c("chr18:67952374-67965471","STAT1")
ProjDir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/"
distFlank<-10000

for(i in 1:length(Genes_names)){
  tryCatch(
    expr = {
      request_test<-atacTracks(distflank = distFlank,Genes = Genes_names[i],Track_list = Track_list,Track_cols = Track_cols,Genome = "hg38")
      pdf(file=paste0(ProjDir,Genes_names[i],".pdf"))
      viewTracks(request_test$tracks, gr=request_test$geneRegion, viewerStyle=request_test$view)
      dev.off()
    },
    error = function(e){ 
      message(paste("No data available for the gene: ",Genes_names[i],sep = ""))
    }
  )
}
