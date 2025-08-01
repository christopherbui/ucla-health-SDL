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

module load cellranger/6.1.1
echo " "

################################################################

CELLRANGER=/u/local/apps/cellranger/6.1.1/cellranger
ID=run_count_1kpbmcs
TRANSCRIPTOME=/u/scratch/c/cbui/cellranger_tutorial/refdata-gex-GRCh38-2020-A
FASTQS=/u/scratch/c/cbui/cellranger_tutorial/pbmc_1k_v3_fastqs
SAMPLE=pbmc_1k_v3

