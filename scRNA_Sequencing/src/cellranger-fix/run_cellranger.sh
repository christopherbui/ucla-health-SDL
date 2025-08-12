#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load perl
module load cellranger/7.1.0

. /u/local/etc/profile.d/sge.sh

ref_dir="/u/scratch/c/cbui/cellranger_tutorial/refdata-gex-GRCh38-2020-A"

# caseID=$(sed -n 1p $myinput)
# fastqPath=$(sed -n 2p $myinput)
# fastqName="${fastqPath##*/}"  ##remove longest before (#) /
# sampID="${fastqName%%_fastqs}"  ##remove longest after (%) _fastqs
out_dir="${caseID}"_"${sampID}"

cellranger count \
    --disable-ui \
    --id=$out_dir \
    --transcriptome=$ref_dir \
    --fastqs=$fastqPath \
    --include-introns=true \
    --jobmode=sge --maxjobs=60 --jobinterval=10000

##--localcores=8 --localmem=64 --maxjobs=64