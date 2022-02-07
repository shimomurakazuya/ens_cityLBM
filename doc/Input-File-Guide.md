# Setting the input json file
Basically, the json file consists of two major categories ```JOB_settings``` and ```citylbm_settings```.  
The former corresponds to the job setting like MPI paralellization, account name and the like.
The latter gives the physical and numerical settings for citylbm.

## ```JOB_settings```
| Variable Name | Type | Explanations | 
| --- | --- | --- |
| `ACCOUNT` | String | Account name or group name. <br> The `ACCOUNT` is often used to charge the requested resources on Supercomputers. |
| `JOB_NAME` | String | Serves as a job name and output directory name | 
| `NB_PROCS` | Int or list of Ints | Number of MPI processes. <br> Equal to number of GPUs if `USE_GPUs` is `true` (flat MPI). |
| `NB_THREADS` | Int | Number of OpenMP threads. Recommended to set as <br> `NB_THREADS` = `NB CPU sockets` x `NB cores per CPU socket` / `NB_PROCS`|
| `OUT_DIR` | String | The output data are stored in this directory (automatically created if not exists). <br> E.g. `/work`, `/data` |
| `TIME` | String | Walltime in hh:mm:ss, mm:ss, or mm. E.g. `qsub -l h_rt=hh:mm:ss` | 
| `USE_GPUs` | Bool | `true` for GPU computing and `false` for CPU computing. |
| `Scaling`  | String |`Strong` or `Weak` (optional). <br> `NB_PROCS` must be given in list of ints with ascending order. See [an example script](https://github.com/hasegawa-yuta-jaea/citylbm/tree/master/run/jobs/example_Tsubame3.0_CPU_strong_scaling.json). |
      
## ```citylbm_settings```
| Variable Name | Type | Explanations | 
| --- | --- | --- |
| `gpu_per_node` | Int | Number of MPI processes per node. E.g. `4` for `SGI8600` supercomputer, `8` for `Aquarius`. |
| `number_of_grid_point` | Int\[2\] | Complicated parameter. Detailed description will be added in [wiki](https://github.com/hasegawa-yuta-jaea/citylbm/wiki/configure_resolution_domain)|
| `domain_min` | Float\[3\] | Origin coordinates (x0, y0, z0) \[meter\] |
| `domain_length` | Float\[3\] | Domain length (Lx, Ly, Lz) \[meter\] |
| `time_end` | Float | Maximum duration of time integral \[minute\] |
| `wrf_start` | Int | Begging time-step of WRF data in `io/input/oklahoma`. <br>`0` = 23:00 CDT, July 16, 2003; 1 step = 1 minute. |
| `velocity_lbm` | Float\[2\] | `[0]` is reference velocity \[m/s\]. <br>`[1]` is ratio of reference velocity and lattice speed.　(lattice speed = Δx/Δt) |
| `cfr_steps` | Int\[2\] | Interval of outputs\[step\]. <br>`[0]` is interval of console output. <br>`[1]` is interval of file output. |
| `restart_flags_and_step` | Int\[2\] | Parameters for restarting calculation from the checkpoint(optional). <br>`[0]` is the flag which enables restarting when `1` is set. <br>`[1]` is the step number where the checkpoint exists. |
| `n_ensemble_members` | Int | Number of ensemble members (optional). |
| `ofs_ensemble_idx` | Int | Offset for the ID of the ensemble members (optional). |
| `cal_rank_div_algo` | String | Algorithm of domain partitioning (optional). <br>Only `"2d_block"` is valid in current.  |
| `cal_rank_div_2d_block` | Int\[2\] | Tiling size for the rank map of 2D domain paritioning (optional). <br> |
