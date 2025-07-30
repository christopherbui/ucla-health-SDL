#!/bin/bash
#$ -cwd
# error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
### change time and resources as needed:
#$ -l h_rt=4:00:00,h_data=64G,exclusive
# Email address to notify (do not modify):
#$ -M $USER@mail
# Notify when
#$ -m bea

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

. /u/local/Modules/default/init/modules.sh
## Specify the cellranger version if you want use a version other than v4.0.0
module load cellranger
module li
echo " "

### change "--help" with your actual cellranger call below:
echo "cellranger --help"
cellranger --help

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "

