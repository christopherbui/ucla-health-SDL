## Example from: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct
## mkdir $SCRATCH/CR_count; cd $SCRATCH/CR_count
## wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
## tar -xvf pbmc_1k_v3_fastqs.tar >> /dev/null 2>&1
## wget https://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
## tar -zxvf refdata-cellranger-GRCh38-3.0.0.tar.gz >> /dev/null 2>&1
## cp /u/local/apps/submit_scripts/submit_cellranger_count_ex_cluster_mode.sh ./ 
## qsub submit_cellranger_count_ex_cluster_mode.sh
#!/bin/bash
#$ -cwd
#  error = Merged with joblog
#$ -o joblog.$JOB_ID
#$ -j y
### change the resources below for the 
### central pipeline as needed (the pipeline
### will dispatch jobs with guessed amount of
### resources independently from what you set
### below):
#$ -l h_data=8G,h_rt=24:00:00
# #$ -pe shared 2
### change the name of the job as/if needed and
### uncomment the line below:
# #$ -N CR_count_ex
#  Email address to notify (please do not change):
#$ -M $USER@mail
#  Notify at beginning and end of job
# #$ -m n   # for never
#$ -m bea

. /u/local/Modules/default/init/modules.sh
### add the specific cellrange version if not using current default version on H2
module load cellranger

. /u/local/etc/profile.d/sge.sh

# echo job info on joblog:
echo "Job $JOB_ID started on:   " `hostname -s`
echo "Job $JOB_ID started on:   " `date `
echo " "

### change [...] with your actual cellranger call below:
echo "cellranger count --id=run_count_1kpbmcs \
--fastqs=$SCRATCH/CR_count/pbmc_1k_v3_fastqs \
--sample=pbmc_1k_v3 \
--transcriptome=$SCRATCH/CR_count/refdata-cellranger-GRCh38-3.0.0 \
--jobmode=sge --maxjobs=60 --jobinterval=10000"
cellranger count --id=run_count_1kpbmcs --fastqs=$SCRATCH/CR_count/pbmc_1k_v3_fastqs --sample=pbmc_1k_v3 --transcriptome=$SCRATCH/CR_count/refdata-cellranger-GRCh38-3.0.0 --jobmode=sge --maxjobs=60 --jobinterval=10000

# echo job info on joblog:
echo "Job $JOB_ID ended on:   " `hostname -s`
echo "Job $JOB_ID ended on:   " `date `
echo " "
