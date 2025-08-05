#!/bin/bash
#$ -cwd
#$ -o /u/scratch/c/cbui/cellranger_tutorial/job-logs/joblog.$JOB_ID
#$ -j y
### CHANGE RESOURCES AS NEEDED:
#$ -l h_rt=3:00:00,h_data=8G,exclusive
#$ -pe shared 4
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

module load cellranger/6.1.1
echo " "

################################################################

ID=/u/scratch/c/cbui/cellranger_tutorial/run_count_1kpbmc
TRANSCRIPTOME=/u/scratch/c/cbui/cellranger_tutorial/refdata-gex-GRCh38-2020-A
FASTQS=/u/scratch/c/cbui/cellranger_tutorial/pbmc_1k_v3_fastqs
SAMPLE=pbmc_1k_v3
JOBMODE=sge

# cellranger count
cellranger count \
    --disable-ui \
    --id=$ID \
    --transcriptome=$TRANSCRIPTOME \
    --fastqs=$FASTQS \
    --sample=$SAMPLE \
    --jobmode=$JOBMODE \
    --jobinterval=10000


################################################################

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "