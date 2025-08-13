#!/bin/bash

. /u/local/Modules/default/init/modules.sh

module load spaceranger/3.1.1

. /u/local/etc/profile.d/sge.sh

ref_dir="/u/scratch/c/cbui/spaceranger_tutorial/all_refs/refdata-gex-mm10-2020-A"

fastqName="${fastqPath##*/}"
sampID="${fastqName%%_*}"
outDir="${caseID}"_"${sampID}"

# no include-introns available
spaceranger count \
    --id="$outDir" \
    --transcriptome="$ref_dir" \
    --probe-set="$csvPath" \
    --fastqs="$fastqPath" \
    --image="$imgPath" \
    --slide="V11J26-127" \
    --area="B1" \
    --create-bam=true \
    --reorient-images=true \
    --disable-ui \
    --jobmode=sge --maxjobs=60 --jobinterval=10000