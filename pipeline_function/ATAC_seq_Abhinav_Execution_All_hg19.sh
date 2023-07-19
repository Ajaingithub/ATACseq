#!/bin/bash

#print the time and date
echo Starting: $(date)

myPROJDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg19"
mySbatchDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg19/sbatch"
mySampleFile="${myPROJDIR}/samplesheet.txt"
myStartModule="all"

##-----------------------------##
##run Sbatch      ##
##-----------------------------##
function runSbatch()
{
SbatchCMD="qsub ${mySbatchDIR}/ATAC_seq_all_Abhinav_JobscriptModular_Longer.sh ${jobSample} ${myStartModule} -o ${jobSample}_log.txt -e ${jobSample}_error.txt"
echo ${SbatchCMD}
eval ${SbatchCMD}
}


#-------------------------------#
#starting job specification     #
#-------------------------------#
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mySampleArray=($(cat ${mySampleFile}))' #sample-index mapping

for i in "${mySampleArray[@]}"
do
jobSample=$i
echo ${jobSample}
runSbatch ${jobSample}
done

#print the time and date again
echo Ending: $(date)
