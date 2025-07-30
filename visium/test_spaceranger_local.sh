#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
### CHANGE RESOURCES AS NEEDED:
#$ -l h_rt=4:00:00,h_data=64G,exclusive
### CHANGE NAME OF JOB AS NEEDED:
## $ -N NAMEOFJOB
### EMAIL ADDRESS TO NOTIFY:
#$ -M $USER@mail
### NOTIFY WHEN
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load environment
. /u/local/Modules/default/init/modules.sh

# load spaceranger
module load spaceranger
echo " "

# spaceranger call
echo "spaceranger --help"
spaceranger --help

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "