#!/bin/bash
#$ -cwd
#  error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
### change the resources below for the 
### central pipeline as needed (the pipeline
### will dispatch jobs with guessed amount of
### resources independently from what you set
### below):
#$ -l h_data=8G,h_rt=24:00:00
# #$ -pe shared 2
### change the name of the job as/if needed and
### uncomment the line below:
# #$ -N NAMEOFYOURJOBIFANY
#  Email address to notify (please do not change):
#$ -M $USER@mail
#  Notify at beginning and end of job
# #$ -m n   # for never
#$ -m bea

. /u/local/Modules/default/init/modules.sh
### add the specific cellrange version if not using current default version on H2
module load cellranger

. /u/local/etc/profile.d/sge.sh

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

### change [...] with your actual cellranger call below:
echo "cellranger [...] --jobmode=sge --maxjobs=60 --jobinterval=10000"
cellranger [...] --jobmode=sge --maxjobs=60 --jobinterval=10000

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
