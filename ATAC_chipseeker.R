+#### Extracting nearby genes with ChipSeeker ######
## To annotate the differential peaks ##
## Run it on the terminal for comparecluster

chipseeker_ATAC <- function(saveDir,dds){
  
}
#setwd("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/pre_EBV_vs_VZV_Old/Differential_all")

# since we ahve already made the differential files during the homer that we will use for the chipseeker


cluster_diff_peaks <- list.files(paste(saveDir,"Homer",sep = ""),pattern = "_coordinates_df.bed",full.names = T)
cluster_diff_peaks_name <- basename(cluster_diff_peaks)
dir.create(paste(saveDir,"peak_annotation",sep = ""),showWarnings = FALSE)
i=0
for (i in 1:length(cluster_diff_peaks)) {
  cluster_peaks <- getPeaks(cluster_diff_peaks[i], sort_peaks=TRUE)
  cluster_peaks_annotated <- annotatePeak(cluster_peaks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, 
                                          tssRegion=c(-1000, 1000), verbose=FALSE)
  
  pdf(paste(saveDir,"peak_annotation/",cluster_diff_peaks_name[i],"_Summary_annotated_peak.pdf",sep = ""), width = 7, height = 4)
  print(plotAnnoPie(cluster_peaks_annotated))
  dev.off()
  
  ## looking at the Enrich pathway using ReactomePA
  pathway_reactome <- enrichPathway(as.data.frame(cluster_peaks_annotated)$geneId, pvalueCutoff=0.05)
  
  pdf(paste(saveDir,"peak_annotation/",cluster_diff_peaks_name[i],"_reactome_pathway.pdf",sep = ""), width = 12, height = 8)
  print(dotplot(pathway_reactome))
  dev.off()
  
  cluster_genes <- bitr(cluster_peaks_annotated@anno$geneId, fromType = "ENTREZID",
                        toType = "SYMBOL", OrgDb = "org.Hs.eg.db")

  write.table(cluster_genes,paste(saveDir,"peak_annotation/",cluster_diff_peaks_name[i],"genes_id_and_name.txt",sep = ""), 
              sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cluster_peaks_annotated@anno$geneName <- cluster_genes[match(cluster_peaks_annotated@anno$geneId, cluster_genes$ENTREZID),2]
  saveRDS(cluster_peaks_annotated, paste(saveDir,"peak_annotation/",cluster_diff_peaks_name[i],".RDS",sep = ""))
}


# Convert the EntrezID to Symbol

#### Performing GO for MF using compareCluster ####
# In order to universal we have to do the Chipseeker for all the peaks in EBV
# at this location /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/pre_Old_vs_young_ebv/Differential
dds <- readRDS(paste(savedir,"dds.RDS",sep = ""))
write.table(rownames(assay(dds)),paste(saveDir,"peak_annotation/all_peaks.txt",sep = ""),
            row.names = FALSE, col.names = FALSE,
            sep = "\t", quote = FALSE)

# sed 's/:/\t/g' all_peaks.txt | sed 's/-/\t/g' | sort -V > all_peaks.bed
all_peaks <- getPeaks(paste(saveDir,"peak_annotation/all_peaks.bed",sep = ""), sort_peaks=TRUE)

all_peaks_annotated <- annotatePeak(all_peaks, TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene, tssRegion=c(-1000, 1000), verbose=FALSE)
all_peaks_annotated_df <- as.data.frame(all_peaks_annotated)

upregulated_geneid <- as.data.frame(unique(up_DF_peaks_annotated@anno$geneId))
upregulated_geneid$Group <- "Upregulated"
colnames(upregulated_geneid) <- c("ENTREZID","Group")

downregulated_geneid <- as.data.frame(unique(down_DF_peaks_annotated@anno$geneId))
downregulated_geneid$Group <- "Downregulated"
colnames(downregulated_geneid) <- c("ENTREZID","Group")
up_and_down_entrez <- rbind(upregulated_geneid,downregulated_geneid)

formula_res_enrichGO_BP_pvalue_cutoff <- compareCluster(ENTREZID~Group, 
                                                        data=up_and_down_entrez, 
                                                        fun="enrichGO",
                                                        OrgDb = "org.Hs.eg.db",
                                                        universe = all_peaks_annotated@anno$geneId,
                                                        ont           = "BP",
                                                        pAdjustMethod = "BH",
                                                        #pvalue = 1,
                                                        pvalueCutoff  = 0.05,
                                                        qvalueCutoff  = 0.05)

formula_res_enrichGO_2 <- formula_res_enrichGO_BP_pvalue_cutoff@compareClusterResult
write.table(formula_res_enrichGO_2, "formula_res_enrichGO_2_pvalue_less_than_0_5.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# formula_res_enrichGO_2_signaling <- formula_res_enrichGO_2[grep("signaling", formula_res_enrichGO_2$Description), ]
# formula_res_enrichGO_2_signaling_Pvalue <- formula_res_enrichGO_2_signaling[formula_res_enrichGO_2_signaling$pvalue < 0.05,]

library("dplyr")
library("tidyr")
formula_res_enrichGO_plot_2 <- separate(formula_res_enrichGO_2, GeneRatio, c("GeneRatio1","GeneRatio2"), sep = "/")
formula_res_enrichGO_plot_2$GeneRatio1<- as.numeric(formula_res_enrichGO_plot_2$GeneRatio1)
formula_res_enrichGO_plot_2$GeneRatio2<- as.numeric(formula_res_enrichGO_plot_2$GeneRatio2)
formula_res_enrichGO_plot_2$GeneRatio <- formula_res_enrichGO_plot_2$GeneRatio1/formula_res_enrichGO_plot_2$GeneRatio2
#formula_res_enrichGO_plot_2 <- formula_res_enrichGO_plot_2[-11,] # since integrin mediated signaling pathway has very low pvalue
formula_res_enrichGO_plot_2$Description <- factor(formula_res_enrichGO_plot_2$Description, 
                                                  levels = formula_res_enrichGO_plot_2$Description[order(formula_res_enrichGO_plot_2$Group)])
library(ggplot2)
p   <-   ggplot(formula_res_enrichGO_plot_2, aes(x=Group,y=Description, color=p.adjust, size=Count)) + 
  geom_point() +
  ggtitle("EBV VS VZV Pvalue < 0.05") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 9),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey")
  )

pdf("EnrichGO_BP_chipseeker_pvalue_05.pdf", width = 8, height = 9)
p
dev.off()

# Manually we have take out the gene from BP for immune, autophagy, cell proliferation, and viral cycle in Google spreadsheet
# ATACSeq Pre and Post vaccination, BP_2KB
# grep -i "upregulated" formula_res_enrichGO_2_pvalue_less_than_0_5.txt | cut -f10 | tr -s "/" "\n" | sort | uniq > upregulated_geneid.txt
# 

# Extracting genes that are upregulated as they have genes which are involved in immune signaling
up_geneid_unique <- read.table("upregulated_geneid.txt")
up_genename_unique <-  up_genes[match(up_geneid_unique$V1,up_genes$ENTREZID, nomatch = 0),2]
write.table(up_genename_unique, 
            "up_genename_unique_EBV.txt", row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")

## Extracting gene for the viral genome inhibition
down_geneid_unique <- read.table("downregulated_geneid.txt")
down_genename_unique <-  down_genes[match(down_geneid_unique$V1,down_genes$ENTREZID, nomatch = 0),2]
write.table(up_genename_unique, 
            "up_genename_unique_2KB_immune.txt", row.names = FALSE,
            col.names = TRUE, quote = FALSE, sep = "\t")

library(session)
save.session("GO_annotation_2KB.Rda")


### Since it worked greatly with p.adjust 0.05. I think we donot need to do the pvalue cut off to be 1
## Now considering withour pvalue cutoff
formula_res_enrichGO_BP <- compareCluster(ENTREZID~Group, 
                                          data=up_and_down_entrez,
                                          fun="enrichGO",
                                          OrgDb = "org.Hs.eg.db",
                                          universe = all_peaks_annotated@anno$geneId,
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          #pvalue = 1,
                                          pvalueCutoff  = 1,
                                          qvalueCutoff  = 1)

formula_res_enrichGO_2 <- formula_res_enrichGO_BP@compareClusterResult
write.table(formula_res_enrichGO_2, "formula_res_enrichGO_2.txt", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)

formula_res_enrichGO_2_signaling <- formula_res_enrichGO_2[grep("signaling", formula_res_enrichGO_2$Description), ]
formula_res_enrichGO_2_signaling_Pvalue <- formula_res_enrichGO_2_signaling[formula_res_enrichGO_2_signaling$pvalue < 0.05,]

library("dplyr")
library("tidyr")
formula_res_enrichGO_plot_2 <- separate(formula_res_enrichGO_2_signaling_Pvalue, GeneRatio, c("GeneRatio1","GeneRatio2"), sep = "/")
formula_res_enrichGO_plot_2$GeneRatio1<- as.numeric(formula_res_enrichGO_plot_2$GeneRatio1)
formula_res_enrichGO_plot_2$GeneRatio2<- as.numeric(formula_res_enrichGO_plot_2$GeneRatio2)
formula_res_enrichGO_plot_2$GeneRatio <- formula_res_enrichGO_plot_2$GeneRatio1/formula_res_enrichGO_plot_2$GeneRatio2
formula_res_enrichGO_plot_2 <- formula_res_enrichGO_plot_2[-23,] # since integrin mediated signaling pathway has very low pvalue
formula_res_enrichGO_plot_2$Description <- factor(formula_res_enrichGO_plot_2$Description, 
                                                  levels = formula_res_enrichGO_plot_2$Description[order(formula_res_enrichGO_plot_2$Group)])
library(ggplot2)
p   <-   ggplot(formula_res_enrichGO_plot_2, aes(x=Group,y=Description, color=p.adjust, size=Count)) + 
  geom_point() +
  ggtitle("EBV VS VZV Pvalue < 0.05") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 9),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey")
  )

pdf("EnrichGO_BP_chipseeker_pvalue_05.pdf", width = 8, height = 9)
p
dev.off()