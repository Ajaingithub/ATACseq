#!/bin/bash
#$ -cwd
#$ -o './'
#$ -e './'
#$ -l h_vmem=8G
#$ -pe threaded 2
#$ -q 1-day
#$ -m bae -M jain.abhinav@mayo.edu

#print the time and date
echo Starting: $(date)

#-------------------------------#
#initialize job context         #
#-------------------------------#
source ~/.bashrc
module load bedtools/2.27.1
module load python/2.7.10

myPROJDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/Ines/CD8_lytic/tohg38/Analysis/downstream/CreateTracks/"
mySampleFile="${myPROJDIR}/SampleFiles.txt"
mySampleFile2="${myPROJDIR}/SampleDepth.txt"

##-----------------------------##
##create FolderLayout          ##
##-----------------------------##
function createFolderLayout()
{
###be aware of the input parameter
###${1} --- means the ${jobSample}

myBedDIR="${myPROJDIR}/bed"
if [ ! -d "${myBedDIR}" ]
then
echo "${myBedDIR} does not exist, creating it now."
tmpCMDSTR="mkdir ${myBedDIR}"
eval "${tmpCMDSTR}"
else
echo "${myBedDIR} already exist, just use it."
fi

myBedgraphDIR="${myPROJDIR}/bedgraph"
if [ ! -d "${myBedgraphDIR}" ]
then
echo "${myBedgraphDIR} does not exist, creating it now."
tmpCMDSTR="mkdir ${myBedgraphDIR}"
eval "${tmpCMDSTR}"
else
echo "${myBedgraphDIR} already exist, just use it."
fi

myBigwigDIR="${myPROJDIR}/bigwig"
if [ ! -d "${myBigwigDIR}" ]
then
echo "${myBigwigDIR} does not exist, creating it now."
tmpCMDSTR="mkdir ${myBigwigDIR}"
eval "${tmpCMDSTR}"
else
echo "${myBigwigDIR} already exist, just use it."
fi
}
##-----------------------------##
##run Bam2Bed                 ##
##-----------------------------##
function runBam2Bed()
{
    tempRead="${jobSample}"
    tempSample="$(basename ${jobSample} | cut -d. -f1)"
    Bam2BedCMD="bamToBed -i ${tempRead} > ${myPROJDIR}/bed/${tempSample}.bed"
    echo ${Bam2BedCMD}
    eval ${Bam2BedCMD}
}

##-----------------------------##
##run Bed truncation to insert ##
##-----------------------------##
function runTruncation()
{
    tempSample="$(basename ${jobSample} | cut -d. -f1)"
    tempRead="${myPROJDIR}/bed/${tempSample}.bed"
    TruncationCMD="awk 'BEGIN{ OFS=\"\\t\"; }{if ("\$6"==\"+\") { start="\$2";end="\$2"+1 } else {start="\$3"-1;end="\$3"}; print "\$1", "start", "end", "\$4", "\$5", "\$6"; }' ${tempRead} | sortBed -i | genomeCoverageBed -bg -i stdin -g ${myPROJDIR}/hg38.chrom.sizes > ${myPROJDIR}/bedgraph/${tempSample}.bedgraph"
    echo ${TruncationCMD}
    eval ${TruncationCMD}
}

##-----------------------------##
##run bedgraph2bigwig          ##
##-----------------------------##
function runNormalization()
{
    tempSample="$(basename ${jobSample} | cut -d. -f1)"
    tempRead="${myPROJDIR}/bedgraph/${tempSample}.bedgraph"
    tempdepth=${depth}
    NormalizationCMD="awk 'BEGIN {OFS=\"\\t\";}{normalized=(("\$4"/${tempdepth})*1000000);print "\$1", "\$2", "\$3", normalized;}' ${tempRead} > ${myPROJDIR}/bedgraph/${tempSample}_normailzed.bedgraph"
    echo ${NormalizationCMD}
    eval ${NormalizationCMD}
}


##-----------------------------##
##run bedgraph2bigwig          ##
##-----------------------------##
function runbedgraph2bigwig()
{
    tempSample="$(basename ${jobSample} | cut -d. -f1)"
    tempRead="${myPROJDIR}/bedgraph/${tempSample}_normailzed.bedgraph"
    bedgraph2bigwigCMD="${myPROJDIR}/bedGraphToBigWig ${tempRead} ${myPROJDIR}/hg38.chrom.sizes ${myPROJDIR}/bigwig/${tempSample}_normailzed.bw"
    echo ${bedgraph2bigwigCMD}
    eval ${bedgraph2bigwigCMD}
}

##-----------------------------##
##run bedgraph2bigwig          ##
##-----------------------------##
function runSmoothing()
{
    tempSample="$(basename ${jobSample} | cut -d. -f1)"
    tempRead="${myPROJDIR}/bigwig/${tempSample}_normailzed.bw"
    SmootingCMD="python ${myPROJDIR}/pyBWsmooth.py -a ${tempRead} -g ${myPROJDIR}/hg38.chrom.sizes"
    echo ${SmootingCMD}
    eval ${SmootingCMD}
}

#-------------------------------#
#starting job specification     #
#-------------------------------#
jobSample=$1
depth=$2
echo "Current Sample is: ${jobSample} - Depth: ${depth}"
echo ${jobSample}
echo ${depth}
createFolderLayout ${jobSample}
runBam2Bed ${jobSample}
runTruncation ${jobSample}
runNormalization ${jobSample}
runbedgraph2bigwig ${jobSample}
runSmoothing ${jobSample}

#print the time and date again
echo Ending: $(date)
