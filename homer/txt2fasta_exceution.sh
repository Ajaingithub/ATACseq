echo starting: $(date)

mySampleFile="samplesheet.txt"
myStartModule="all"

## Run sbatch ##
function runSbatch()
{
SbatchCMD="qsub txt2fasta_bash.sh ${jobSample} ${myStartModule} -o ${jobSample}_log.txt -e ${jobSample}_error.txt"
echo ${SbatchCMD}
eval ${SbatchCMD}
}

## starting the job specification
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mySampleArray=($(cat ${mySampleFile}))' #sample-index mapping

for i in "${mySampleArray[@]}"
do
jobSample=$i
echo ${jobSample}
runSbatch ${jobSample}
done
