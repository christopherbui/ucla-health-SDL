#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load perl
module load cellranger/7.1.0

. /u/local/etc/profile.d/sge.sh

# change reference path as necessary
ref_dir="/u/scratch/c/cbui/cellranger_tutorial/refdata-gex-GRCh38-2020-A"

# results_dir="/u/scratch/c/cbui/cellranger_tutorial/results"
# mkdir -p "${results_dir}"
# cd "${results_dir}"

# caseID=$(sed -n 1p $myinput)
# fastqPath=$(sed -n 2p $myinput)
fastqName="${fastqPath##*/}"
sampID="${fastqName%%_*}"
outDir="${caseID}"_"${sampID}"

cellranger count \
    --id=$outDir \
    --transcriptome=$ref_dir \
    --fastqs=$fastqPath \
    --include-introns=true \
    --disable-ui \
    --jobmode=sge --maxjobs=60 --jobinterval=10000

##--localcores=8 --localmem=64 --maxjobs=64