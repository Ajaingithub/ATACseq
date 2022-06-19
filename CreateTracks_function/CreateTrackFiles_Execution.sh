#!/bin/bash

#print the time and date
echo Starting: $(date)

myPROJDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/CreateTracks"
mySbatchDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/Analysis/Downstream_Run7_Good_QC_filter_removed_2_Outlier_adding_3_doubtful_add_24_samples_include_VZV/sbatch"
mySampleFile="${myPROJDIR}/SampleFiles.txt"
mySampleFile2="${myPROJDIR}/SampleDepth.txt"

##-----------------------------##
##run Sbatch      ##
##-----------------------------##
function runSbatch()
{
SbatchCMD="qsub ${mySbatchDIR}/CreateTrackFiles.sh ${jobSample} ${depth} -o ${jobSample}_log.txt -e ${jobSample}_error.txt"
echo ${SbatchCMD}
eval ${SbatchCMD}
}


#-------------------------------#
#starting job specification     #
#-------------------------------#
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mySampleArray=($(cat ${mySampleFile}))' #sample-index mapping
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mySampleArray2=($(cat ${mySampleFile2}))' #sample-index mapping 
index=0

for i in "${mySampleArray[@]}"
do
jobSample=$i
depth=${mySampleArray2[$index]}
echo ${jobSample}
echo ${depth}
runSbatch ${jobSample} ${depth}
index=$index+1
done

#print the time and date again
echo Ending: $(date)
