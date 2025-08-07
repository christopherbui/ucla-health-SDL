# Hoffman2 Computing

`qrsh -l h_rt=4:00:00,h_data=4G`

- `qrsh` requests resources to start interactive session (processes run on interactive node have access to requested resources)

`qsub script.sh`

- `qsub` submits a job to the scheduler

You can specify resources as an argument for `qsub`, or specify it in the script.

```bash
qsub -l h_rt=200,h_data=100M -o joblog -j y add_two_vectors_sequentially.sh
```



`qstat` - check job status

```sh
[cbui@login2 ~]$ qstat -u $USER
job-ID     prior   name       user         state submit/start at     queue             jclass                   slots ja-task-ID
-----------------------------------------------------------------------------------------------------------
  10052492 0.00000 test_space cbui         qw    07/30/2025 15:00:13                                                           2
```

`qstat -j <jobID>` - see any error messages for why job failed

```sh
[cbui@login1 job-logs]$ qstat -j 10094350
==============================================================
job_number:                 10094350
jclass:                     NONE
exec_file:                  job_scripts/10094350
submission_time:            08/05/2025 14:58:34.295
owner:                      cbui
uid:                        22520
group:                      sdubinet
gid:                        10646
supplementary group:        sdubinet
sge_o_home:                 /u/home/c/cbui
sge_o_log_name:             cbui
sge_o_path:                 /u/local/apps/cellranger/6.1.1/cellranger-6.1.1/bin:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/bin/intel64:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/bin:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/bin:/u/local/compilers/intel/2020.4/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin:/u/local/compilers/intel/2020.4/debugger_2020/gdb/intel64/bin:/u/systems/UGE8.6.4/bin/lx-amd64:/u/local/bin:/u/local/sbin:/u/local/Modules/4.7.0/gcc-4.8.5/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/u/home/c/cbui/bin
sge_o_shell:                /bin/bash
sge_o_workdir:              /u/scratch/c/cbui/cellranger_tutorial/script
sge_o_host:                 login1
account:                    sge
cwd:                        /u/scratch/c/cbui/cellranger_tutorial/script
merge:                      y
hard resource_list:         exclusive=TRUE,h_data=8G,h_rt=10800
mail_options:               abe
mail_list:                  cbui@mail
notify:                     FALSE
job_name:                   submit_cellranger_cluster.sh
stdout_path_list:           NONE:NONE:/u/scratch/c/cbui/cellranger_tutorial/job-logs/joblog.$JOB_ID
priority:                   0
jobshare:                   0
restart:                    n
env_list:                   SCRATCH=/u/scratch/c/cbui
script_file:                submit_cellranger_cluster.sh
parallel environment:       shared range: 4
project:                    sdubinet_prj
department:                 defaultdepartment
binding:                    NONE
mbind:                      NONE
submit_cmd:                 qsub submit_cellranger_cluster.sh
category_id:                372
request_dispatch_info:      FALSE
start_time            1:    12/31/1969 16:00:00.000
job_state             1:    Eqw
error reason    1:          08/05/2025 18:12:26 [22520:15235]: execvp(/work/UGE/8.6.4/n7679/job_scripts/10094350, "/work/UGE/8.6.4/n7679/job_scripts/10094350") failed: No such file or directory
scheduling info:            (Collecting of scheduler job information is turned off)

```



`id`

`myresources`

`jobs`

`qacct -j $jobID`

- shows data for job



## Requesting Resources

```bash
# start interactive session
NUMcores=5; GBperCORE=4; TIME="04:00:00"; qrsh -pe shared "${NUMcores}" \
-l i,h_rt="${TIME}",h_data="${GBperCORE}",h_vmem="$((NUMcores*GBperCORE))"G
```



## Multiple Jobs

Each `qsub` submission specifies its resources.

Use a wrapper script that dynamically:

- calls `qsub` with specific resources
- applies string manipulation to files & outputs

*SCRIPT.sh* takes in 2 positional arguments.

```bash
#!/bin/bash
TASKnCores="${1}"; inFQ="${2}"
OUTbase="$(basename "${inFQ}")"
OUTbase="${OUTbase%.fastq.gz}".BSBolt-GRCh38Ens114

module load stuff

# do stuff; ${1} & ${2} used here
```



```bash
#!/bin/bash
TASKnCores=2; TASKperCoreGB=7; TASKmaxTIME="1:00:00"
for x in SRS*-head-part??.fastq.gz; do	# SRX101207se-head-part01.fastq.gz
	SMPshrt="${x%%-*}"			# SRX101207se
	SMPshrt="${SMPshrt#SRX}"	# 101207se
	SMPshrt="${SMPshrt%se}"		# 101207
	PARTnum="${x%.fastq.gz}"	# SRX101207se-head-part01
	PARTnum="${PARTnum##*-part}"		# 01
	JOBname="A${SMPshrt}_${PARTnum}"	# A101207_01
	
	qsub -cwd -r no -j no -m as -N "${JOBname}" -pe shared "${TASKnCores}" \
	-l "h_rt=${TASKmaxTIME},h_data=${TASKperCoreGB}G,h_vmem=$((TASKnCores*TASKperCOREGB))G" \
	SCRIPT.sh "${TASKnCores}" "${x}"
done
```



## Estimate Resource Usage

C

## Other Shell

`awk`

`pv`

`cut`

`gzip -c -d`

`paste - - - -`



# Test Job Submission

`qacct -j 10110504`

```bash
[cbui@login3 pbmc_1k_v3_fastqs]$ qacct -j 10110504
==============================================================
qname        pod_smp.q
hostname     n1064
group        sdubinet
owner        cbui
project      sdubinet_prj
department   defaultdepartment
jobname      submit_cellranger_cluster_LF.sh
jobnumber    10110504
taskid       undefined
pe_taskid    NONE
account      sge
priority     0
cwd          /u/scratch/c/cbui/cellranger_tutorial/script
submit_host  login3
submit_cmd   qsub submit_cellranger_cluster_LF.sh
qsub_time    08/07/2025 13:20:43.535
start_time   08/07/2025 13:24:51.439
end_time     08/07/2025 14:59:18.705
granted_pe   shared
slots        4
failed       0
deleted_by   NONE
exit_status  0
ru_wallclock 5667.266
ru_utime     10132.795
ru_stime     838.089
ru_maxrss    15698988
ru_ixrss     0
ru_ismrss    0
ru_idrss     0
ru_isrss     0
ru_minflt    328648496
ru_majflt    1269
ru_nswap     0
ru_inblock   97338880
ru_oublock   52296368
ru_msgsnd    0
ru_msgrcv    0
ru_nsignals  0
ru_nvcsw     40042567
ru_nivcsw    58901189
wallclock    5667.561
cpu          10970.884
mem          131510.327
io           416.474
iow          99.270
ioops        145029429
maxvmem      17.697G
maxrss       0.000
maxpss       0.000
arid         undefined
jc_name      NONE
bound_cores  NONE

```



