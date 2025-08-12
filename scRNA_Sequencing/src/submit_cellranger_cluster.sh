#!/bin/bash

myinput="/u/scratch/c/cbui/cellranger_tutorial/script/cr_input.txt"


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


# submit qsub cellranger job
qsub \
    -cwd \
    -o "/u/scratch/c/cbui/cellranger_tutorial/job-logs" \
    -j y \
    -l h_rt=24:00:00,h_data=8G \
    -pe shared 8 \
    -M $USER@mail \
    -m bea \
    -v myinput=$myinput \
    run_cellranger.sh