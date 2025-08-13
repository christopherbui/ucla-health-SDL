#!/bin/bash

# change as necessary
myinput="/u/scratch/c/cbui/cellranger_tutorial/script/myinput.txt"

# change as necessary
results_dir="/u/scratch/c/cbui/cellranger_tutorial/results"
mkdir -p "$results_dir"


###################################################################
# check input file existence
if [[ ! -f "$myinput" ]]; then
    echo "ERROR: No Valid Input File: "
    echo "$myinput"
    exit 1
fi

# check input file contents
echo "--------------------------------"
cat "$myinput"
echo "--------------------------------"

# ask for user confirmation
read -p "Continue with job submission? (y/n): " answer
case "$answer" in
    [yY]* )
        echo "Submitting job..."
        ;;
    [nN]* )
        echo "Cancelling..."
        exit 0
    ;;
    * )
        echo "Invalid input. Use (y/n)..."
        exit 1
        ;;
esac
###################################################################

# for each project name & fastq location, read 2 lines
total_lines=$(wc -l < "${myinput}")
if (( total_lines %2 != 0 )); then
    echo "ERROR: Input File Has Uneven Number of Lines..."
    exit 1
fi

while read caseID && read fastqPath; do
    fastqName="${fastqPath##*/}"
    sampID="${fastqName%%_*}"
    outDir="${caseID}"_"${sampID}"

    # check if fastq directory exists
    if [[ ! -d "$fastqPath" ]]; then
        echo "ERROR: FASTQ Directory Not Found..."
        echo "$fastqPath"
        exit 1
    fi

    echo "caseID: $caseID"
    echo "fastqPath: $fastqPath"
    echo "fastqName: $fastqName"
    echo "sampID: $sampID"
    echo "outDir: $outDir"


    # submit qsub cellranger job
    qsub \
        -wd "$results_dir" \
        -o "/u/scratch/c/cbui/cellranger_tutorial/job-logs" \
        -j y \
        -l h_rt=24:00:00,h_data=8G \
        -pe shared 8 \
        -M $USER@mail \
        -m bea \
        -v caseID="$caseID",fastqPath="$fastqPath" \
        run_cellranger.sh

done < "$myinput"