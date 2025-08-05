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

