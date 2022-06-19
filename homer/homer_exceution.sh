echo starting: $(date)

mySampleFile="homer_samplesheet.txt"
myBackgroundFile="background_samplesheet.txt"

### Run Sbatch #####
function runSbatch()
{
SbatchCMD="qsub homer_bash.sh ${jobSample} ${backgroundSample} -o ${jobSample}_${backgroundSample}_log.txt -e ${jobSample}_${backgroundSample}_error.txt"
echo ${SbatchCMD}
eval ${SbatchCMD}
}

## starting the job specification
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mySampleArray=($(cat ${mySampleFile}))' #sample-index mapping
IFS=$'\r\n' GLOBIGNORE='*' command eval  'mybackgroundArray=($(cat ${myBackgroundFile}))' #background-index mapping

paste -d@ homer_samplesheet.txt background_samplesheet.txt | while IFS="@" read -r f1 f2
do
jobSample="$f1"
backgroundSample="$f2"
echo ${jobSample}
echo ${backgroundSample}
runSbatch ${jobSample} ${backgroundSample}
done


