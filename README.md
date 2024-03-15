# ATACseq
Automation of the ATACseq pipeline

## Preprocessing

The bash script start from fastqs to Peaks calling. Generated for hg19 and hg38 genome build. The preprocessing pipeline is developed in collaboration with Dr. Rohit R. Jadhav

#### The process include
1. Fastq QC
2. Trimming
3. Alignment
4. Alignment QC
5. Peak Calling


        sh pipeline_function/ATAC_seq_Abhinav_Execution_All_hg19.sh # Invoke pipeline_function/ATAC_seq_all_Abhinav_JobscriptModular_Longer_hg19.sh
        sh pipeline_function/ATAC_seq_Abhinav_Execution_All_hg38.sh # Invoke pipeline_function/ATAC_seq_all_Abhinav_JobscriptModular_Longer_hg38.sh

This pipeline will submit multiple samples into the cluster compute nodes and will run parallely. 

## Postprocessing
This include filtering of the ATAC peaks and samples to be used for downstream analysis
1. Make the Splot that could be used for identify the overallping peaks and merging of the samples

         source("./pipeline_function/splot_ATACseq.R")
         splot(peak_path=sampleinfo2$peak_location,
               saveDir="./output/",
               build = "hg38",
               species = "human")

2. Based on the ouput from the splot_ATACseq.R. Filter out the non-overlap peaksm blacklisted region peaks, and low quality peaks so that the matrix could be further use for the downstream analysis

         source("./pipeline_function/After_Splot.R")
         After_Splot(All_merge_path = All_merge_path,
                     sampleinfo =  sampleinfo2,
                     overlaps = 5,
                     saveDir = savedir,
                     bam_path = sampleinfo2$bam_location,
                     blacklist_region = hg38_blacklist_region)

## Downstream Analysis
This analysis will perform PCA, Differential, kmeans clustering, ChromVar, Homer, and bigwig merging.
1. **PCA**
   
        source("./pipeline_function/PCA.R")
        PCA_plot(cnts_path = cnts_path,
                 sampleinfo = sampleinfo,
                 saveDir = savedir,
                 group_level = group_level,
                 label = "Person",
                 color = "Group",
                 shape = "Batch")
   
2. **Differential**: Based on the PCA result perform differential
   
       source("./pipeline_function/differential_ATAC.R")
       ATAC_differential(mainDir = maindir,
                         saveDir = savedir,
                         compare_group = "Group",
                         comparison=comparison,
                         main_name = "CD8_EBV")

3. **Kmeans clustering & heatmap**:
Using the differential result we develop the kmeans clustering
a. **kmeans_clustering_ATAC** : This function will generate the gap statistics to identify how many cluster (K) can be generated

           source("pipeline_function/k_mer_clustering.R")
           kmeans_clustering_ATAC(differential_peaks = diff_peaks_vec,
                                   rlog_dds = rld,
                                   sampleinfo = sampleinfo,
                                   Group = "Group",
                                   saveDir = savedir)
   
b. **k_means_plots.R** : This function will generate the heatmap with foreground and background peaks for each cluster so that we can use HOMER to identify the Enriched Transcription factor activity within the cluster.

        source("pipeline_function/k_means_plots.R")
        k_means_plot(cluster_gap = gap_stat_2, 
                        df_scale = df_scale, 
                        saveDir = savedir, 
                        clus_num = 14,
                        deseq_dataset = deseq_dataset)

4. **Homer** :
   We will generate the Enriched Transcription factor output from each of the clusters for which we have foreground and background peaks
   
           sh homer/homer_exceution.sh # Invoke this script homer/homer_bash.sh

   This pipeline will submit multiple Homer for each kmeans cluster into the cluster compute nodes and will run parallely.

   Now we will merge K_means_plot with the Homer Output

           source("./pipeline_function/Heatmap_and_Cluster_Homer.R")
           homer_ATAC(saveDir = savedir,
                      final_2 = final_2,
                      df_scale = df_scale,
                      colindex_arrange = colindex_arrange,
                      column_anno = column_anno)

   The Output can be visualize here ATACseq/Rplots.pdf

6. Then further we performed the ChipSeeker and rGREAT to annotate the peaks
           pipeline_function/ATAC_Great.R

Also you can perform preprocessing step using these function

The Version Control has been done through RStudio

If there is any issue please comment.
