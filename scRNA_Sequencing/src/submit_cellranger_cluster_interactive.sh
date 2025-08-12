#!/bin/bash

echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

# load environment
. /u/local/Modules/default/init/modules.sh

module load cellranger/6.1.1
echo " "

################################################################

ID=run_count_1kpbmc
TRANSCRIPTOME=/u/scratch/c/cbui/cellranger_tutorial/refdata-gex-GRCh38-2020-A
FASTQS=/u/scratch/c/cbui/cellranger_tutorial/pbmc_1k_v3_fastqs
SAMPLE=pbmc_1k_v3
JOBMODE=local

# cellranger count
cellranger count \
    --disable-ui \
    --id=$ID \
    --transcriptome=$TRANSCRIPTOME \
    --fastqs=$FASTQS \
    --sample=$SAMPLE \
    --jobmode=$JOBMODE \
    --jobinterval=10000 \
    --localcores=4 \
    --localmem=32


################################################################

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
