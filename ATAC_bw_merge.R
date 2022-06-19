+
  ### Merging bw files 
  ### Merging bigwig files for ebv pre old cm, em, 
  ##### BW and merge the BW ###########
### Most of the stuff needs to be done in Bash but we can make the bash script using R commmands
## Ran this on the bash scripting making the bed --> bedgraph --> bigwig files.
## the code /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/sbatch 
## from the location you can identify the data. Keep h_vmem = 10G, and -pe threaded 10
## qsub -pe threaded 10 -l h_vmem=10G -N ATAC_seq -q 1-hour -o ${mySbatchDIR}/${tempSample}_log.txt -e ${mySbatchDIR}/${tempSample}_error.txt ${mySbatchDIR}/CreateTrackFiles.sh ${jobSample} ${depth}

bw_merge <- function(trackDir,saveDir,){
  
}
trackDir <- "/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg38/Analysis/downstream/CreateTracks/bigwig/bigwig_files/"
setwd(trackDir)
bwFileName <- list.files(pattern = ".bw.s20.w150sw.bw", full.names = TRUE)
bwFilepath <- list.files(path = trackDir,pattern = ".bw.s20.w150sw.bw", full.names = TRUE)
sampleinfo <- read.table("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg38/Analysis/downstream/sampleinfo2.txt",
                         header = T, sep = "\t")

source("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Resources/R_functions_Abhinav/partial_matching_Abhinav.R")
sampleinfo$bwFilepath <- bwFilepath[partial_match_Abhinav(sampleinfo$SampleID2,bwFilepath)]
dplyr::select(sampleinfo,SampleID2,bwFileName,Group)

write.table(sampleinfo,paste(saveDir,"sampleinfo3.txt",sep = ""),
            sep = "\t",quote = F, row.names = F, col.names = T)

# for comparison 1 and comparison 2 we are merging everything
rm(first_line)
rm(merge_cmd)
first_line <- vector()
merge_cmd <- vector()
groups <- unique(sampleinfo[,Group])
for(i in 1:length(groups)){
  print(grep(groups[i],sampleinfo[,2]))
  bw_file_name <- sampleinfo[grep(groups[i],sampleinfo[,Group]),'bwFileName']
  first_line[i] <- paste("mySampleGroup=\"",groups[i],"\"", sep = "")
  merge_cmd[i] <- paste("mergeCMD=\"qsub -pe threaded 10 -l h_vmem=10G -N ATAC_seq -q 1-hour -o ./log_files/${mySampleGroup}_log.txt -e ./error_files/${mySampleGroup}_error.txt -cwd ${myScriptDIR}/merge_bigwig.sh -g 10 -T ${myFileDIR} ${myFileDIR}/${mySampleGroup}_Average.bw ${myPROJDIR}/hg38.chrom.sizes",
                        bw_file_name[1],bw_file_name[2],bw_file_name[3],bw_file_name[4], bw_file_name[5], bw_file_name[6], bw_file_name[7], 
                        bw_file_name[8],bw_file_name[9],bw_file_name[10],bw_file_name[11],bw_file_name[12],bw_file_name[13], bw_file_name[14], 
                        bw_file_name[15], bw_file_name[16],bw_file_name[17],bw_file_name[18],bw_file_name[19],bw_file_name[20],bw_file_name[21], 
                        bw_file_name[22], bw_file_name[23],bw_file_name[24],bw_file_name[25],bw_file_name[26],bw_file_name[27],bw_file_name[28],
                        bw_file_name[29], bw_file_name[30], 
                        "\"", sep = " ")
}

final_merge <- vector()
for(i in 1:length(merge_cmd)){
  final_merge[i] <- paste(first_line[i],merge_cmd[i],"echo ${mergeCMD}", "eval ${mergeCMD}","" , sep = "\n")
}

cat(final_merge)
setwd("/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg38/Analysis/downstream/CreateTracks/bigwig")
write.table(final_merge,"final_merge.txt",sep=" ", quote = FALSE, row.names = FALSE, col.names = FALSE)
system("sed -i 's: NA::g' final_merge.txt")
system("sed -i 's:w150sw.bw \":w150sw.bw\":g' final_merge.txt" )
system ("cat directory_head.txt final_merge.txt > EBV_bw_merge.sh")
system("sh EBV_bw_merge.sh")














