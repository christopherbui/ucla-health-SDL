#!/bin/bash

. /u/local/Modules/default/init/modules.sh

module load spaceranger/3.1.1

. /u/local/etc/profile.d/sge.sh

ref_dir="/u/scratch/c/cbui/spaceranger_tutorial/refdata-gex-mm10-2020-A"
img_dir="/u/scratch/c/cbui/spaceranger_tutorial/datasets/Visium_FFPE_Mouse_Brain_image.jpg"

out_dir="${caseID}"_"${sampID}"

# no include-introns available
spaceranger count \
    --id=$out_dir \                                
    --transcriptome=$ref_dir \
    --create-bam=true \
    --image=$img_dir \
    --disable-ui