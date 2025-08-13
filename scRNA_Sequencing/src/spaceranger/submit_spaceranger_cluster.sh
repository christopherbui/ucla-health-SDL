#!/bin/bash

# change as necessary
myinput="/u/scratch/c/cbui/spaceranger_tutorial/script/myinput.txt"

# change as necessary
results_dir="/u/scratch/c/cbui/spaceranger_tutorial/results"
mkdir -p "$results_dir"

# datasets directory with image & probe set .csv
datasets_dir="/u/scratch/c/cbui/spaceranger_tutorial/datasets"


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

# CHECK: for each project name & fastq location, read 2 lines
total_lines=$(wc -l < "${myinput}")
if (( total_lines %2 != 0 )); then
    echo "ERROR: Input File Has Uneven Number of Lines..."
    exit 1
fi

while read caseID && read fastqPath; do
    # check if fastq directory exists
    if [[ ! -d "$fastqPath" ]]; then
        echo "ERROR: FASTQ Directory Not Found..."
        echo "$fastqPath"
        exit 1
    fi

    fastqName="${fastqPath##*/}"
    sampID="${fastqName%%_*}"
    outDir="${caseID}"_"${sampID}"

    # get image
    imgPath=$(find "${datasets_dir}" -maxdepth 1 -type f -name "${sampID}_*.jpg")
    # get probe set csv
    csvPath=$(find "${datasets_dir}" -maxdepth 1 -type f -name "${sampID}_*.csv")
    # CHECK: should only be 1 image and 1 csv for each sampID
    imgCount=$(echo "${imgPath}" | grep -c .)
    csvCount=$(echo "${csvPath}" | grep -c .)
    if (( imgCount != 1 )); then
        echo "ERROR: Expected exactly 1 image..."
        echo "${imgPath}"
        exit 1
    fi
    if (( csvCount != 1 )); then
        echo "ERROR: Expected exactly 1 csv..."
        echo "${csvPath}"
        exit 1
    fi        

    echo "caseID: $caseID"
    echo "fastqPath: $fastqPath"
    echo "fastqName: $fastqName"
    echo "sampID: $sampID"
    echo "outDir: $outDir"

    # submit qsub spaceranger job
    qsub \
        -wd "$results_dir" \
        -o "/u/scratch/c/cbui/spaceranger_tutorial/job-logs" \
        -j y \
        -l h_rt=24:00:00,h_data=8G \
        -pe shared 8 \
        -M $USER@mail \
        -m bea \
        -v caseID="$caseID",fastqPath="$fastqPath",imgPath="$imgPath",csvPath="$csvPath" \
        run_spaceranger.sh

done < "$myinput"