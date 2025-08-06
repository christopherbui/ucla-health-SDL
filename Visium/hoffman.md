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

