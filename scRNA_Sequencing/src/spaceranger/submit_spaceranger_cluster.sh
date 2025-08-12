#!/bin/bash

myinput="/u/scratch/c/cbui/spaceranger_tutorial/myinput.txt"


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

for (( i=1; i <= total_lines; i+=2 )); do
    caseID=$(sed -n "${i}p" "${myinput}")
    fastqPath=$(sed -n "$((i+1))p" "${myinput}")
    fastqName="${fastqPath##*/}"
    sampID="${fastqName%%_fastqs}"
    out_dir="${caseID}"_"${sampID}"

    # check if fastq directory exists
    if [[ ! -d "$fastqPath" ]]; then
        echo "ERROR: FASTQ Directory Not Found..."
        echo "$fastqPath"
    fi

    echo "Lines ${i} & $((i+1)):"
    echo "caseID: $caseID"
    echo "fastqPath: $fastqPath"
    echo "fastqName: $fastqName"
    echo "sampID: $sampID"
    echo "out_dir: $out_dir"
    echo ""

    # submit qsub cellranger job
    qsub \
        -cwd \
        -o "/u/scratch/c/cbui/spaceranger_tutorial/job-logs" \
        -j y \
        -l h_rt=24:00:00,h_data=8G \
        -pe shared 8 \
        -N ${caseID}_${sampID} \
        -M $USER@mail \
        -m bea \
        -v caseID=$caseID,fastqPath=$fastqPath,sampID=$sampID \
        run_spaceranger.sh
done