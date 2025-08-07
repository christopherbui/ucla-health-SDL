#!/bin/bash
## cellranger count (10X genome) for 3' primer platform
## input is text file list (1) project name (e.g. CIRM, GGO, CCL21) at the first line,
## and (2) location of fastq files at the rest. All fastq files with
## same pre-fix (i.e. before underscore) of sampleID in the folder 
## will be analyzed together by cellranger count.
## output is in location that job is submitted. Name: projectName_sampleID
## requires 8GB x8 shared for at least 24hrs.
## Check readme file for re-run as job killed before finishing  
## Use cellranger v7.1.0 with default mode including intron and cluster mode
## 8/31/23

. /u/local/Modules/default/init/modules.sh
module load perl
module load cellranger/7.1.0

. /u/local/etc/profile.d/sge.sh

ref_dir="/u/project/sdubinet/linhtran/refdata-gex-mm10-2020-A"

caseID=$(sed -n 1p $myinput)
fastqF=$(sed -n 2p $myinput)
fastqName="${fastqF##*/}"  ##remove longest before (#) /
fastqPath="${fastqF%/*}"
sampID="${fastqName%%_*}"  ##remove longest after (%) _
out_dir="${caseID}"_"${sampID}"

cellranger count --id=$out_dir \
--transcriptome=$ref_dir \
--fastqs=$fastqPath \
--include-introns=true \
--jobmode=sge --maxjobs=60 --jobinterval=10000

##--localcores=8 --localmem=64 --maxjobs=64

