#!/usr/bin/env bash
# alignAll.sh
# Usage: bash scripts/alignAll.sh 1>results/logs/alignAll.log 2>results/logs/alignAll.err &

# Initialize variable to contain the directory of un-trimmed fastq files 
fastqPath="/work/courses/BINF6309/AiptasiaMiSeq/fastq/"

# Initialize variable to contain the suffix for the left reads
leftSuffix=".R1.fastq"
rightSuffix=".R2.fastq"

outDir='quant/'

function alignAll {
    # Loop through all the left-read fastq files in $fastqPath
    for leftInFile in $fastqPath*$leftSuffix
    do
    # Remove the path from the filename and assign to pathRemoved
    pathRemoved="${leftInFile/$fastqPath/}"
    # Remove the left-read suffix from $pathRemoved and assign to suffixRemoved
    sample="${pathRemoved/$leftSuffix/}"
    # Print $sampleName to see what it contains after removing the path
    echo sample
    salmon quant -l IU \
        -1 $fastqPath${sample}.R1.fastq \
        -2 $fastqPath${sample}.R2.fastq \
        -i AipIndex \
        --validateMappings \
        -o ${outDir}${sample}
    done
}
alignAll 
