#!/bin/bash
#$ -cwd
#$ -o './'
#$ -l h_vmem=16G
#$ -pe threaded 10
#$ -q 1-day
#$ -m bae -M jain.abhinav@mayo.edu

echo Starting: $(date)
source ~/.bashrc

function running_homer()
{
echo "Running Homer"
echo "${jobSample}"
basename "${jobSample}"

# Getting the basename of the file
f="$(basename -- $jobSample)"
echo "$f"

# Removing the extension from the file
g="${f%.*}"
echo "$g"

homer_cmd="findMotifs.pl ${g}.fa human ./${g}/ -fasta ${background}"
echo "${homer_cmd}"
eval "${homer_cmd}"
}

# Running jobs
jobSample="${1}"
background="${2}"
running_homer ${jobSample}

echo Ending: $(date)
