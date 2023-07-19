#!/bin/bash
#$ -cwd
#$ -o '/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38/sbatch/qsub_output'
#$ -l h_vmem=16G
#$ -pe threaded 2
#$ -q 1-day
#$ -m bae -M jain.abhinav@mayo.edu

#print the time and date
echo Starting: $(date)

#-------------------------------#
#initialize job context         #
#-------------------------------#
source ~/.bashrc
# module load python/3.9.2
module load fastqc/0.11.8
#module load bowtie2/2.3.3.1
module load picard/2.21.6
module load java/14.0.1 samtools/1.10 bedtools/2.27.1 r/R-4.0.3
# module load legacy ngsplot/2.47
module load python/2.7.10

myDATADIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/raw_data"
myPROJDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/toHg38"
mySCRIPTSDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/ATAC_scripts"
mypreseqDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/preseq_v2.0"
myigvDIR="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/ATAC_seq/IGV_2.10.3"

##-----------------------------##
##create FolderLayout          ##
##-----------------------------##
function createFolderLayout()
{
###be aware of the input parameter
###${1} --- means the ${jobSample}

myAnalyDIR="${myPROJDIR}/Analysis"
if [ ! -d "${myAnalyDIR}" ]
then
echo "${myAnalyDIR} does not exist, creating it now."
tmpCMDSTR="mkdir ${myAnalyDIR}"
eval "${tmpCMDSTR}"
else
echo "${myAnalyDIR} already exist, just use it."
fi


mytrimDIR="${myAnalyDIR}/trim"
if [ ! -d "${mytrimDIR}" ]
then
echo "${mytrimDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${mytrimDIR}; mkdir ${mytrimDIR}/${1}"
tmpCMDSTR="mkdir ${mytrimDIR}; mkdir ${mytrimDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${mytrimDIR} already exist, just use it."
#####for each sample#####
#sampletrimDir=${mytrimDIR}/${1}
sampletrimDir=${mytrimDIR}/${jobSample}
if [ ! -d "${sampletrimDir}" ]
then
echo "${sampletrimDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${sampletrimDir}"
eval "${tmpCMDSTR}"
else
echo "${sampletrimDir} already exist, just use it."
fi
fi

myfastqcDIR="${myAnalyDIR}/fastqc"
if [ ! -d "${myfastqcDIR}" ]
then
echo "${myfastqcDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${myfastqcDIR}; mkdir ${myfastqcDIR}/${1}"
tmpCMDSTR="mkdir ${myfastqcDIR}; mkdir ${myfastqcDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${myfastqcDIR} already exist, just use it."
#####for each sample#####
#samplefastqcDir=${myfastqcDIR}/${1}
samplefastqcDir=${myfastqcDIR}/${jobSample}
if [ ! -d "${samplefastqcDir}" ]
then
echo "${samplefastqcDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${samplefastqcDir}"
eval "${tmpCMDSTR}"
else
echo "${samplefastqcDir} already exist, just use it."
fi
fi

mybowtieDIR="${myAnalyDIR}/bowtie"
if [ ! -d "${mybowtieDIR}" ]
then
echo "${mybowtieDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${mybowtieDIR}; mkdir ${mybowtieDIR}/${1}"
tmpCMDSTR="mkdir ${mybowtieDIR}; mkdir ${mybowtieDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${mybowtieDIR} already exist, just use it."
#####for each sample#####
#samplebowtieDir=${mybowtieDIR}/${1}
samplebowtieDir=${mybowtieDIR}/${jobSample}
if [ ! -d "${samplebowtieDir}" ]
then
echo "${samplebowtieDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${samplebowtieDir}"
eval "${tmpCMDSTR}"
else
echo "${samplebowtieDir} already exist, just use it."
fi
fi

myQCDIR="${myAnalyDIR}/QC"
if [ ! -d "${myQCDIR}" ]
then
echo "${myQCDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${myQCDIR}; mkdir ${myQCDIR}/${1}"
tmpCMDSTR="mkdir ${myQCDIR}; mkdir ${myQCDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${myQCDIR} already exist, just use it."
#####for each sample#####
#sampleQCDir=${myQCDIR}/${1}
sampleQCDir=${myQCDIR}/${jobSample}
if [ ! -d "${sampleQCDir}" ]
then
echo "${sampleQCDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${sampleQCDir}"
eval "${tmpCMDSTR}"
else
echo "${sampleQCDir} already exist, just use it."
fi
fi

myPeaksDIR="${myAnalyDIR}/Peaks"
if [ ! -d "${myPeaksDIR}" ]
then
echo "${myPeaksDIR} does not exist, creating it now."
#tmpCMDSTR="mkdir ${myPeaksDIR}; mkdir ${myPeaksDIR}/${1}"
tmpCMDSTR="mkdir ${myPeaksDIR}; mkdir ${myPeaksDIR}/${jobSample}"
eval "${tmpCMDSTR}"
else
echo "${myPeaksDIR} already exist, just use it."
#####for each sample#####
#samplePeaksDir=${myPeaksDIR}/${1}
samplePeaksDir=${myPeaksDIR}/${jobSample}
if [ ! -d "${samplePeaksDir}" ]
then
echo "${samplePeaksDir} does not exist, creating it now."
tmpCMDSTR="mkdir ${samplePeaksDir}"
eval "${tmpCMDSTR}"
else
echo "${samplePeaksDir} already exist, just use it."
fi
fi

}

##-----------------------------##
##run Trimming                 ##
##-----------------------------##
function runTrim()
{

tempRead1="${myDATADIR}/${jobSample}*_1.fq"
tempRead2="${myDATADIR}/${jobSample}*_2.fq"

trimCMD="python ${mySCRIPTSDIR}/pyadapter_trim.py -a ${tempRead1} -b ${tempRead2}"
echo ${trimCMD}
eval ${trimCMD}
}

##-----------------------------##
##move Trim Files              ##
##-----------------------------##
function mvTrimFile()
{
myTrimDIR="${myPROJDIR}/sbatch"

tempTrimFile1="${myTrimDIR}/${jobSample}*_1.trim.fastq"
tempTrimFile2="${myTrimDIR}/${jobSample}*_2.trim.fastq"

mvTrimFileCMD1="mv ${tempTrimFile1} ${myPROJDIR}/Analysis/trim/${jobSample}/"
echo ${mvTrimFileCMD1}
eval ${mvTrimFileCMD1}

mvTrimFileCMD2="mv ${tempTrimFile2} ${myPROJDIR}/Analysis/trim/${jobSample}/"
echo ${mvTrimFileCMD2}
eval ${mvTrimFileCMD2}

}

##-----------------------------##
##run Fastqc                   ##
##-----------------------------##
function runfastqc()
{

tempRead1="${myDATADIR}/${jobSample}*_1.fq"
tempRead2="${myDATADIR}/${jobSample}*_2.fq"

tempTrim1="${myPROJDIR}/Analysis/trim/${jobSample}/${jobSample}*_1.trim.fastq"
tempTrim2="${myPROJDIR}/Analysis/trim/${jobSample}/${jobSample}*_2.trim.fastq"
fastqcCMD="fastqc -o ${myPROJDIR}/Analysis/fastqc/${jobSample} -f fastq ${tempRead1} ${tempRead2}"
echo ${fastqcCMD}
eval ${fastqcCMD}
fastqcCMD_Trim="fastqc -o ${myPROJDIR}/Analysis/fastqc/${jobSample}_trim -f fastq ${tempTrim1} ${tempTrim2}"
echo ${fastqcCMD_Trim}
eval ${fastqcCMD_Trim}
}

##-----------------------------##
##run Bowtie                   ##
##-----------------------------##
function runbowtie()
{
myTrimDIR="${myPROJDIR}/Analysis/trim/${jobSample}"

tempRead1="${myTrimDIR}/${jobSample}*_1.trim.fastq"
tempRead2="${myTrimDIR}/${jobSample}*_2.trim.fastq"

bowtieCMD="/research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/miniconda3/bin/bowtie2 --no-discordant --no-dovetail --no-mixed --no-unal --very-sensitive -X 2000 -p 25 -x /research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/RefGenomes/H_sapiens/hg38/indexes/bowtie2/hg38 --rg-id ${jobSample} -1 ${tempRead1} -2 ${tempRead2} -S ${myPROJDIR}/Analysis/bowtie/${jobSample}/${jobSample}.sam"
echo ${bowtieCMD}
eval ${bowtieCMD}
}

##-----------------------------##
##run QualityControl           ##
##-----------------------------##
function runQC()
{
myBowtieDir="${myPROJDIR}/Analysis/bowtie/${jobSample}"
mybaseQCDIR="${myPROJDIR}/Analysis/QC"
myQCDIR="${myPROJDIR}/Analysis/QC/${jobSample}"

#converts a .sam to .bam file (more compressed) along with removing reads with mapping quality < 20
samtobamCMD="samtools view -q 20 -uT /research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/RefGenomes/H_sapiens/hg38/genome.fa ${myBowtieDir}/${jobSample}.sam -o ${myQCDIR}/${jobSample}.bam"
echo ${samtobamCMD}
eval ${samtobamCMD}

#sorts file to be in order; samples need to be sorted before they can be merged
bamsortCMD="samtools sort -@ 6  -T ${myQCDIR}/tmp  -o ${myQCDIR}/${jobSample}.sorted.bam ${myQCDIR}/${jobSample}.bam"
echo ${bamsortCMD}
eval ${bamsortCMD}

#generates an index file for later processing step
bamindexCMD="samtools index ${myQCDIR}/${jobSample}.sorted.bam"
echo ${bamindexCMD}
eval ${bamindexCMD}

#saves reads from chrom 1-22 in a file (Remove chrM, chrX, chrY)
rmchrCMD="samtools view -u ${myQCDIR}/${jobSample}.sorted.bam `for c in $(seq 1 22); do echo -n "chr${c} "; done` | samtools rmdup - ${myQCDIR}/${jobSample}.rmdup.bam"
echo ${rmchrCMD}
eval ${rmchrCMD}

#removes duplicate reads
#rmdupCMD="samtools rmdup ${myQCDIR}/${jobSample}.sorted.bam ${myQCDIR}/${jobSample}.rmdup.bam"
#echo ${rmdupCMD}
#eval ${rmdupCMD}

#Sort reads again
sortrmdup="samtools sort -@ 6 -T ${myQCDIR}/tmp -o ${myQCDIR}/${jobSample}.rmdup.sorted.bam ${myQCDIR}/${jobSample}.rmdup.bam"
echo ${sortrmdup}
eval ${sortrmdup}

#generates index for later processing steps
rmdupindexCMD="samtools index ${myQCDIR}/${jobSample}.rmdup.sorted.bam"
echo ${rmdupindexCMD}
eval ${rmdupindexCMD}

#converts the bam file to a bedgraph file
bamtobedgraphCMD="genomeCoverageBed -bg -split -ibam ${myQCDIR}/${jobSample}.rmdup.sorted.bam -trackopts 'type=bedGraph name='${jobSample} -g /research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/RefGenomes/H_sapiens/hg38/genome.fa.fai > ${myQCDIR}/${jobSample}.rmdup.bedGraph"
echo ${bamtobedgraphCMD}
eval ${bamtobedgraphCMD}

# convert to tdf
bamtotdf="${myigvDIR}/igvtools count ${myQCDIR}/${jobSample}.rmdup.sorted.bam ${myQCDIR}/${jobSample}.rmdup.tdf /research/labs/immunology/goronzy_weyand/GoronzyLab_SCG/RefGenomes/H_sapiens/hg38/genome.fa.fai"
echo ${bamtotdf}
eval ${bamtotdf}

#QC metrics
#Maps TSS openness as a measurement of atac-seq quality
#vplotCMD="python ${mySCRIPTSDIR}/pyMakeVplot.py -a ${myQCDIR}/${jobSample}.rmdup.sorted.bam -b ${mySCRIPTSDIR}/hglft_parsed_hg38_RefSeq.merged.bed -p ends -e 2000 -u -v -c 6 -o ${myQCDIR}/${jobSample}.bed.vect"
#echo ${vplotCMD}
#eval ${vplotCMD}

#generates fragement plot, also a measurement of openness
#insertsizeCMD="picard CollectInsertSizeMetrics INPUT=${myQCDIR}/${jobSample}.rmdup.sorted.bam OUTPUT=${myQCDIR}/${jobSample}.InsertSizeMetrics.txt HISTOGRAM_FILE=${myQCDIR}/${jobSample}.InsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT"
#echo ${insertsizeCMD}
#eval ${insertsizeCMD}

#estimates library complexity
#preseqcurveCMD="${mypreseqDIR}/preseq c_curve -o ${myQCDIR}/${jobSample}.c.curve.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
#echo ${preseqcurveCMD}
#eval ${preseqcurveCMD}
#preseqextrapCMD="${mypreseqDIR}/preseq lc_extrap -o ${myQCDIR}/${jobSample}.lc.extrap.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
#echo ${preseqextrapCMD}
#eval ${preseqextrapCMD}
}

##-----------------------------##
##run QC Metrics Seprately     ##
##-----------------------------##
function runQCMetrics()
{
myBowtieDir="${myPROJDIR}/Analysis/bowtie/${jobSample}"
mybaseQCDIR="${myPROJDIR}/Analysis/QC"
myQCDIR="${myPROJDIR}/Analysis/QC/${jobSample}"

#QC metrics
#Maps TSS openness as a measurement of atac-seq quality
vplotCMD="python ${mySCRIPTSDIR}/pyMakeVplot.py -a ${myQCDIR}/${jobSample}.rmdup.sorted.bam -b ${mySCRIPTSDIR}/hglft_parsed_hg38_RefSeq.merged.bed -p ends -e 2000 -u -v -c 6 -o ${myQCDIR}/${jobSample}.bed.vect"
echo ${vplotCMD}
eval ${vplotCMD}

#generates fragement plot, also a measurement of openness
insertsizeCMD="picard CollectInsertSizeMetrics INPUT=${myQCDIR}/${jobSample}.rmdup.sorted.bam OUTPUT=${myQCDIR}/${jobSample}.InsertSizeMetrics.txt HISTOGRAM_FILE=${myQCDIR}/${jobSample}.InsertSizeMetrics.pdf VALIDATION_STRINGENCY=SILENT"
echo ${insertsizeCMD}
eval ${insertsizeCMD}

#estimates library complexity
preseqcurveCMD="${mypreseqDIR}/preseq c_curve -o ${myQCDIR}/${jobSample}.c.curve.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
echo ${preseqcurveCMD}
eval ${preseqcurveCMD}
preseqextrapCMD="${mypreseqDIR}/preseq lc_extrap -o ${myQCDIR}/${jobSample}.lc.extrap.txt -B ${myQCDIR}/${jobSample}.sorted.bam"
echo ${preseqextrapCMD}
eval ${preseqextrapCMD}
}

##-----------------------------##
##run V-Plot Seprately     ##
##-----------------------------##
function runVPlot()
{
myBowtieDir="${myPROJDIR}/Analysis/bowtie/${jobSample}"
mybaseQCDIR="${myPROJDIR}/Analysis/QC"
myQCDIR="${myPROJDIR}/Analysis/QC/${jobSample}"

#QC metrics
#Maps TSS openness as a measurement of atac-seq quality
vplotCMD="python ${mySCRIPTSDIR}/pyMakeVplot.py -a ${myQCDIR}/${jobSample}.rmdup.sorted.bam -b ${mySCRIPTSDIR}/hglft_parsed_hg38_RefSeq.merged.bed -p ends -e 2000 -u -v -c 6 -o ${myQCDIR}/${jobSample}.bed.vect"
echo ${vplotCMD}
eval ${vplotCMD}
}

##-----------------------------##
##run MACS (Peak calling)      ##
##-----------------------------##
function runMACS()
{
myPeakDIR="${myPROJDIR}/Analysis/Peaks/${jobSample}"

MACSCMD="macs2 callpeak -t ${myPROJDIR}/Analysis/QC/${jobSample}/${jobSample}.rmdup.sorted.bam -f BAM --keep-dup all -n ${jobSample} --nomodel --shift -100 --extsize 200 -q .01 --nolambda --call-summits --outdir ${myPeakDIR}"
echo ${MACSCMD}
eval ${MACSCMD}
}

##-----------------------------##
##run Gunzip      ##
##-----------------------------##
function runGunzip()
{
GZIPCMD1="gunzip ${myDATADIR}/${jobSample}*_1.fq.gz"
GZIPCMD2="gunzip ${myDATADIR}/${jobSample}*_2.fq.gz"
echo ${GZIPCMD1}
eval ${GZIPCMD1}
echo ${GZIPCMD2}
eval ${GZIPCMD2}
}

##-----------------------------##
##run Gunzip      ##
##-----------------------------##
function runGzip()
{
GZIPCMD1="gzip ${myDATADIR}/${jobSample}*_1.fq"
GZIPCMD2="gzip ${myDATADIR}/${jobSample}*_2.fq"
echo ${GZIPCMD1}
eval ${GZIPCMD1}
echo ${GZIPCMD2}
eval ${GZIPCMD2}
}


#-------------------------------#
#starting job specification     #
#-------------------------------#

jobSample=$1
startModule=$2
echo "Current Sample is: ${jobSample}"

if [ ${startModule} = "all" ]
then
    echo "Running all modules.."
    createFolderLayout ${jobSample}
    runGunzip ${jobSample}
    runTrim ${jobSample}
    mvTrimFile ${jobSample}
    runfastqc ${jobSample}
    runbowtie ${jobSample}
    runQC ${jobSample}
    runMACS ${jobSample}
    runQCMetrics ${jobSample}
    runGzip ${jobSample}
elif [ ${startModule} = "bowtie" ]
then
    echo "Running bowtie.."
    runbowtie ${jobSample}
    runQC ${jobSample}
    runMACS ${jobSample}
    runQCMetrics ${jobSample}
    runGzip ${jobSample}
elif [ ${startModule} = "qc" ]
then
    echo "Running QC.."
    runQC ${jobSample}
    runMACS ${jobSample}
    runQCMetrics ${jobSample}
    runGzip ${jobSample}
elif [ ${startModule} = "metrics" ]
then
    echo "Running metrics.."
    runQCMetrics ${jobSample}
    runGzip ${jobSample}
elif [ ${startModule} = "peakonly" ]
then
    echo "Running Peak calling only!!"
    runMACS ${jobSample}
elif [ ${startModule} = "zip" ]
then
    echo "Running GZip!!"
    runGzip ${jobSample}
else
    echo "Please provide appropriate arguments!!"
fi


#print the time and date again
echo Ending: $(date)
