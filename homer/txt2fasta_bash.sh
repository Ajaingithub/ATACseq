#!/bin/bash
#$ -cwd
#$ -o './'
#$ -l h_vmem=16G
#$ -pe threaded 4
#$ -q 1-day
#$ -m bae -M jain.abhinav@mayo.edu

echo Starting: $(date)
source ~/.bashrc

function text2bed() 
{
echo "Running text2bed"
echo "${jobSample}"
basename "${jobSample}"

# Getting the basename of the file
f="$(basename -- $jobSample)"
echo "$f"

# Removing the extension from the file
g="${f%.*}"
echo "$g"

while IFS= read -r line
do
  sed 's/:/\t/g' "$line" | sed 's/-/\t/g' > "$g".bed
done <<<"${jobSample}"
}

function bed2fasta()
{
echo "Running bed2fasta"
echo "${jobSample}"
basename "${jobSample}"

# Getting the basename of the file
f="$(basename -- $jobSample)"
echo "$f"

# Removing the extension from the file
g="${f%.*}"
echo "$g"

bed2fasta_cmd="bedtools getfasta -fo ${g}.fa -fi /research/labs/immunology/goronzy_weyand/GoronzyLab_Mayo/Abhinav/igenome_base/references/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa -bed ${g}.bed"
echo "${bed2fasta_cmd}"
eval "${bed2fasta_cmd}"
}

# Running jobs
jobSample="${1}"
text2bed ${jobSample}
bed2fasta ${jobSample}

echo Ending: $(date)

