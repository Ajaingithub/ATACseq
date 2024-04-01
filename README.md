# ATACseq
Automation of the ATACseq pipeline

The information about each parameter has been provided at the top in the R file of each function

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

![Screenshot 2024-04-01 at 5 58 02 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/7314925b-be5e-4a3a-bd79-4ce8f2175233) ![Screenshot 2024-04-01 at 5 58 28 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/9e3a586d-37f4-427d-b2a0-834f25884819)


## Postprocessing
This include filtering of the ATAC peaks and samples to be used for downstream analysis
1. Make the Splot that could be used for identify the overallping peaks and merging of the samples

         source("./pipeline_function/splot_ATACseq.R")
         splot(peak_path=sampleinfo2$peak_location,
               saveDir="./output/",
               build = "hg38",
               species = "human")

   ![Screenshot 2024-04-01 at 6 01 19 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/45298048-8cbb-4b08-8a60-3fd54766f3de)


3. Based on the ouput from the splot_ATACseq.R. Filter out the non-overlap peaksm blacklisted region peaks, and low quality peaks so that the matrix could be further use for the downstream analysis

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
   
   ![Screenshot 2024-04-01 at 5 59 55 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/d5c2a99f-bd69-4c82-90d8-6c78e3b5b93f)

3. **Differential**: Based on the PCA result perform differential
   
       source("./pipeline_function/differential_ATAC.R")
       ATAC_differential(mainDir = maindir,
                         saveDir = savedir,
                         compare_group = "Group",
                         comparison=comparison,
                         main_name = "CD8_EBV")

4. **Kmeans clustering & heatmap**:
Using the differential result we develop the kmeans clustering
a. **kmeans_clustering_ATAC** : This function will generate the gap statistics to identify how many cluster (K) can be generated

           source("pipeline_function/k_mer_clustering.R")
           kmeans_clustering_ATAC(differential_peaks = diff_peaks_vec,
                                   rlog_dds = rld,
                                   sampleinfo = sampleinfo,
                                   Group = "Group",
                                   saveDir = savedir)
   
   ![Screenshot 2024-04-01 at 6 02 45 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/99c55e7c-3eae-4a40-b7ba-ef70838a2940)

b. **k_means_plots.R** : This function will generate the heatmap with foreground and background peaks for each cluster so that we can use HOMER to identify the Enriched Transcription factor activity within the cluster.

        source("pipeline_function/k_means_plots.R")
        k_means_plot(cluster_gap = gap_stat_2, 
                        df_scale = df_scale, 
                        saveDir = savedir, 
                        clus_num = 14,
                        deseq_dataset = deseq_dataset)

![Screenshot 2024-04-01 at 6 05 37 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/cce8e81b-dd94-4ab7-9195-ff1729e21472)

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
   
![Screenshot 2024-04-01 at 6 04 00 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/c23af4fd-c6c0-44bd-a08a-8315ea9c341f)

6. Then further we performed the ChipSeeker and rGREAT to annotate the peaks

    `pipeline_function/ATAC_Great.R`

8. **ChormVar** : ChromVar also perform the Transcription factor enrichment including all the peaks irrespective it is differential or not. 

        source("./pipeline_function/ATAC_ChromVar.R")

        ATAC_ChromVar(mainDir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg38/Analysis/downstream/",
                      saveDir = "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg38/Analysis/downstream/ChromVar/",
                      sampleinfo = read.table(paste(mainDir,"sampleinfo3.txt", sep = ""),sep = "\t", header = T) # this provide the metadata information,
                      SampleID = SampleId,
                      SampleDepth = "SampleDepth",
                      SampleID_2 = SampleID_2,
                      core=32,
                      min_peaks=0.15,
                      min_depth=1500,
                      Group = "Group",
                      col_arrange = c("naïve","cm","em","temra","lat2","lyt1","lyt2"))
   
![Screenshot 2024-04-01 at 6 06 55 PM](https://github.com/Ajaingithub/ATACseq/assets/37553954/8ec809d6-ffc4-4d39-b712-955883ac5e91)

9. **BigWig Merge** : Combine the bigwig files based on the groups like naive, cm, em etc.

  ` pipeline_function/ATAC_bw_merge.R    `


   
If there is any issue please comment.
