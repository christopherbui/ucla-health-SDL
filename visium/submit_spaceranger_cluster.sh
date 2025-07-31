#!/bin/bash
#$ -cwd
#$ -o joblog.$JOB_ID
#$ -j y
### CHANGE RESOURCES AS NEEDED:
#$ -l h_rt=4:00:00,h_data=16G,exclusive
#$ -pe shared 2
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

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "
# --------------------------------------------------------------


# spaceranger call

# define paths
SPACERANGER=/u/local/apps/spaceranger/2.1.0/spaceranger
TRANSCRIPTOME=refdata-gex-mm10-2020-A
PROBE_SET=Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv
FASTQS=datasets/Visium-FFPE_Mouse_Brain_fastqs
IMAGE=datasets/Visium-FFPE_Mouse_Brain_image.jpg

echo "spaceranger [...] --jobmode=sge --maxjobs=60 --jobinterval=10000"

# Run spaceranger
$SPACERANGER count \
  --id=sample123 \
  --transcriptome=$TRANSCRIPTOME \
  --probe-set=$PROBE_SET \
  --fastqs=$FASTQS \
  --image=$IMAGE \
  --slide=V19J01-123 \
  --area=B1 \
  --jobmode=sge \
  --maxjobs=60 \
  --localcores=2 \
  --localmem=32 \
  --jobinterval=10000



# --------------------------------------------------------------
# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "