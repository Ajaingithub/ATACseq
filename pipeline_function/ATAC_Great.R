### Running GREAT for the annotation fo the peaks in each clysters
# Using rGREAT package for 1st Annotation #########
setwd(paste(saveDir,"Homer",sep = ""))
files <- list.files(paste(saveDir,"Homer",sep = ""),pattern = "_coordinates_df.bed",full.names = T)
files <- files[grep("all", files, invert = TRUE)]
# awk '{print $0"\t"$1":"$2"-"$3}' cnts_noblacklist_data.bed > cnts_noblacklist_data_2.bed
bground <- read.table(paste(saveDir,"cnts_noblacklist_2.bed",sep = ""), sep = "\t")

i=0
for (i in 1:length(files)) {
  bed_file <- read.table(files[i],sep = "\t")
  bed_file$V4 <- paste("C",i,":",bed_file$V1,":",bed_file$V2,"-",bed_file$V3,sep = "")
  filename <- paste("C",i,"_bed_file",sep = "")
  assign(filename,bed_file)
}

bed_file <- ls(pattern = "_bed_file")
jobname <- paste("C",1:length(files),"_GREAT_job",sep = "")

i=0
library("rGREAT")
for (i in 1:length(files)){
  bed <-  get(bed_file[i])
  job = submitGreatJob(bed, bg = bground,
                       species = "hg38",
                       includeCuratedRegDoms = TRUE,
                       rule = "basalPlusExt",
                       adv_upstream = 5.0,
                       adv_downstream = 1.0,
                       adv_span = 1000.0)
  assign(jobname[i], job)
}

### Extracting the GO and genes term for each Clusters
jobs <- ls(pattern="GREAT_job")
jobs <- jobs[grep("job$",jobs)]
for (i in 1:length(jobs)){
  GO_table_name <- paste(jobs[i], "GO_ontology", sep = "_")
  GO_table <- getEnrichmentTables(get(jobs[i]), category = "GO")
  assign(GO_table_name, GO_table)
  gene_table_name <- paste(jobs[i], "Genes_ontology", sep = "_")
  gene_table <- getEnrichmentTables(get(jobs[i]), category = "Genes")
  assign(gene_table_name, gene_table)
}


GO_ontology <- ls(pattern = "_GREAT_job_GO_ontology")

#Molecular Function
for (i in 1:length(GO_ontology)){
  GOMF <- get(GO_ontology[i])$`GO Molecular Function`
  GOMF_selected <- dplyr::select(GOMF, name, Hyper_Adjp_BH, Hyper_Foreground_Region_Hits, Hyper_Background_Gene_Hits)
  colnames(GOMF_selected) <- c("Biological Process", "AdjPvalue", "RegionHits","GeneHits")
  GOMF_selected_top_20 <- head(GOMF_selected,20)
  GOMF_name <- paste(GO_ontology[i],"MF",sep = "_")
  assign(GOMF_name, GOMF_selected_top_20)
}

# Adding the cluster name for each df
MF_files = ls(pattern = "_GREAT_job_GO_ontology_MF")
clusters <- paste("C",1:length(MF_files), sep = "")
for(i in 1:length(MF_files)){
  assign(MF_files[i], `[[<-`(get(MF_files[i]), "Cluster", value = clusters[i]))
}

# Combining all the MF datasets for all the clusters
GOMF_combined = rbind(C1_GREAT_job_GO_ontology_MF, C2_GREAT_job_GO_ontology_MF, C3_GREAT_job_GO_ontology_MF,
                      C4_GREAT_job_GO_ontology_MF, C5_GREAT_job_GO_ontology_MF,C6_GREAT_job_GO_ontology_MF, C7_GREAT_job_GO_ontology_MF) 
# C6_GREAT_job_GO_ontology_MF,
# C7_GREAT_job_GO_ontology_MF, C8_GREAT_job_GO_ontology_MF, C9_GREAT_job_GO_ontology_MF)

GOMF_combined_pvalue <- GOMF_combined[GOMF_combined$AdjPvalue<0.01,]
colnames(GOMF_combined_pvalue) <- c("Molecular_Function", colnames(GOMF_combined_pvalue)[-1])


GOMF_combined_pvalue <- remove.factors(GOMF_combined_pvalue)
GOMF_combined_pvalue$Molecular_Function <- factor(GOMF_combined_pvalue$Molecular_Function, 
                                                  levels = unique(GOMF_combined_pvalue$Molecular_Function))


p <- ggplot(GOMF_combined_pvalue, aes(y=Molecular_Function, x=Cluster, color = AdjPvalue, size = GeneHits)) + 
  geom_point() +
  theme_bw()

dir.create(paste(saveDir,"GO_ontology",sep = ""),showWarnings = FALSE)
setwd(paste(saveDir,"GO_ontology",sep = ""))
pdf("Cluster_molecular_function_2.pdf", width = 12, height = 12)
p
dev.off()

## Biological Process
for (i in 1:length(GO_ontology)){
  GOBP <- get(GO_ontology[i])$`GO Biological Process`
  GOBP_selected <- dplyr::select(GOBP, name, Hyper_Adjp_BH, Hyper_Foreground_Region_Hits, Hyper_Background_Gene_Hits)
  colnames(GOBP_selected) <- c("Biological Process", "AdjPvalue", "RegionHits","GeneHits")
  GOBP_selected_top_20 <- head(GOBP_selected,20)
  GOBP_name <- paste(GO_ontology[i],"BP",sep = "_")
  assign(GOBP_name, GOBP_selected_top_20)
}

# Adding the cluster name for each df
BP_files = ls(pattern = "_GREAT_job_GO_ontology_BP")
BP_files <- BP_files[grep("_GREAT_job_GO_ontology_BP$", BP_files)]
clusters <- paste("C",1:7, sep = "")
for(i in 1:c(1:7)){
  assign(BP_files[i], `[[<-`(get(BP_files[i]), "Cluster", value = clusters[i]))
}

# Combining all the BP datasets for all the clusters
GOBP_combined = rbind(C1_GREAT_job_GO_ontology_BP, C2_GREAT_job_GO_ontology_BP, C3_GREAT_job_GO_ontology_BP,
                      C4_GREAT_job_GO_ontology_BP, C5_GREAT_job_GO_ontology_BP, C6_GREAT_job_GO_ontology_BP, C7_GREAT_job_GO_ontology_BP)
# C6_GREAT_job_GO_ontology_BP,
# C7_GREAT_job_GO_ontology_BP, C8_GREAT_job_GO_ontology_BP, C9_GREAT_job_GO_ontology_BP)

GOBP_combined_pvalue <- GOBP_combined[GOBP_combined$AdjPvalue<0.01,]
colnames(GOBP_combined_pvalue) <- c("Biological_Process", colnames(GOBP_combined_pvalue)[-1])


GOBP_combined_pvalue <- remove.factors(GOBP_combined_pvalue)
GOBP_combined_pvalue$Biological_Process <- factor(GOBP_combined_pvalue$Biological_Process, 
                                                  levels = unique(GOBP_combined_pvalue$Biological_Process))

p <- ggplot(GOBP_combined_pvalue, aes(y=Biological_Process, x=Cluster, color = AdjPvalue, size = GeneHits)) + 
  geom_point() +
  theme_bw()

pdf("Cluster_biological_process.pdf", width = 12, height = 15)
p
dev.off()

### In Biological Process only considering Signaling pathways and pvalue less than 0.05
for (i in 1:length(GO_ontology)){
  GOBP <- get(GO_ontology[i])$`GO Biological Process`
  GOBP_selected <- dplyr::select(GOBP, name, Hyper_Adjp_BH, Hyper_Foreground_Region_Hits, Hyper_Background_Gene_Hits)
  colnames(GOBP_selected) <- c("Biological_Process", "AdjPvalue", "RegionHits","GeneHits")
  GOBP_selected_signaling_Pvalue <- GOBP_selected[grep("signal",GOBP_selected[GOBP_selected$AdjPvalue < 0.05,'Biological_Process']),]
  GOBP_selected_signaling_Pvalue_top20 <- GOBP_selected_signaling_Pvalue[order(-GOBP_selected_signaling_Pvalue$GeneHits),] %>% head(20)
  GOBP_name <- paste(GO_ontology[i],"BP_signaling_significant_top20",sep = "_")
  assign(GOBP_name, GOBP_selected_signaling_Pvalue_top20)
}

BP_files = ls(pattern = "GO_ontology_BP_signaling_significant_top20")
clusters <- paste("C",1:7, sep = "")
for(i in c(1:5,7)){
  assign(BP_files[i], `[[<-`(get(BP_files[i]), "Cluster", value = clusters[i]))
}

GOBP_combined_significant_top20 = rbind(C1_GREAT_job_GO_ontology_BP_signaling_significant_top20, C2_GREAT_job_GO_ontology_BP_signaling_significant_top20, C3_GREAT_job_GO_ontology_BP_signaling_significant_top20,
                                        C4_GREAT_job_GO_ontology_BP_signaling_significant_top20, C5_GREAT_job_GO_ontology_BP_signaling_significant_top20, C7_GREAT_job_GO_ontology_BP_signaling_significant_top20) 
#C8_GREAT_job_GO_ontology_BP_signaling_significant_top20, C9_GREAT_job_GO_ontology_BP_signaling_significant_top20)


GOBP_combined_significant_top20 <- GOBP_combined_significant_top20[GOBP_combined_significant_top20$AdjPvalue < 0.01,]
GOBP_combined_significant_top20 <- remove.factors(GOBP_combined_significant_top20)
GOBP_combined_significant_top20$Biological_Process <- factor(GOBP_combined_significant_top20$Biological_Process, 
                                                             levels = unique(GOBP_combined_significant_top20$Biological_Process))

p <- ggplot(GOBP_combined_significant_top20, aes(y=Biological_Process, x=Cluster, color = AdjPvalue, size = GeneHits)) + 
  geom_point() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y=element_text(size=9),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.minor = element_line(colour = "grey"),
        panel.grid.major = element_line(colour = "grey"))

pdf("Cluster_biological_process_signaling_top_20.pdf", width = 8, height = 8)
p
dev.off()
