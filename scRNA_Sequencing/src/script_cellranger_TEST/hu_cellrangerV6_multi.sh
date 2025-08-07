#!/bin/bash
## cellranger multi (v6.0.2) for 5' and V(D)J together
## --> input includes (1) configuration csv file and (2)
## (2) input for bash script with line 1= SampleID(output folder)
## line 2 =name of csv file
## output is in location where job is submitted
## Note: cellranger options can be specified by configuration csv
## 07/12/21

. /u/local/Modules/default/init/modules.sh
module load perl
export PATH=/u/project/sdubinet/linhtran/cellranger-6.0.2:$PATH

. /u/local/etc/profile.d/sge.sh

###ref_dir="/u/project/sdubinet/linhtran/ref4cellranger_GRch38v300"

sampleID=$(sed -n 1p $myinput)
configFile=$(sed -n 2p $myinput)

cellranger multi --id=$sampleID \
--csv=$configFile \
--jobmode=sge --maxjobs=60 --jobinterval=10000
##--localcores=8 --localmem=64 --maxjobs=60 --jobinterval=10000
##--jobmode=sge --maxjobs=60 --jobinterval=10000
