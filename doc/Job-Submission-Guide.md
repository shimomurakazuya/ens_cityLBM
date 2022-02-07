# Design of job submission helper
Job submission helper is designed to ease job sumissions on different supercomputer platforms.
This document describes the basic usage and porting to new supercomputer platforms.

## Requirements
The job submission helper is written in python. It relys on python 3.4 or later, but does not depend on any external python libraries.

## Usage
Firstly, the environmental variables ```SUPERCOMPUTER``` and ```DEVICE``` should be set in ```~/.bash_profile```. In the following, 
it is assumed that ```SUPERCOMPUTER``` is set as ```<New_supercomputer_name>``` and ```DEVICE``` is set as ```<Device_name>```.
Then, it is needed to place <your_input_file>.json on run/jobs. For example, the job submission on SGI8660 can be made by
```
cd run
python sub.py example_sgi8600.json --dir ./foo
```
The helper automatically generates the job script on the specified environment followed by the job submission. If you place your ```<your_input_file>.json``` on ```run/foo```, the optional argument ```--dir``` (default: ```./jobs```) should be set as ```./foo```.

## Available devices and Supercomputers
Currently, the follwoing devices and supercomputers are supported.
| `DEVICE` | Remarks |
| --- | --- |
| `p100` | Nvidia [Pascal generation GPU](https://images.nvidia.com/content/pdf/tesla/whitepaper/pascal-architecture-whitepaper.pdf). Available on Tsubame3.0. |
| `v100` | Nvidia [Volta generation GPU](https://images.nvidia.com/content/volta-architecture/pdf/volta-architecture-whitepaper.pdf). Available on ABCI, SGI8600, DGX-2. |
| `A100` | Nvidia [Ampere generation GPU](https://images.nvidia.com/aem-dam/en-zz/Solutions/data-center/nvidia-ampere-architecture-whitepaper.pdf). Available on Aquarius. |
| `gnu` | Compiling for conventional CPUs with g++. Basically used for test. |
| `bdw` | Intel Broadwell CPU. Intel compilers icpc and mpiicpc are necessary. Available on Tsubame3.0. |
| `a64fx` | Fujitsu ARM [A64FX CPU](https://github.com/fujitsu/A64FX). Available on Fugaku |
| `a64fx_flow` | Fujitsu ARM [A64FX CPU](https://github.com/fujitsu/A64FX). Available on Flow. |

| `SUPERCOMPUTERS` | Remarks |
| --- | --- |
| `SGI8600` | Cascadelake and V100 machine in JAEA @ Japan |
| `Tsubame3.0` | Broadwell and P100 machine in Tokyo Tech @ Japan |
| `Aquarius` | Icelake and A100 machine in Univ. Tokyo @ Japan |
| `Flow` | A64FX machine in Nagoya Univ. @ Japan |

To run on a new supercomputer ```<New_supercomputer_name>``` equipped with ```<Device_name>```, 
it is need to update the helper in the following manner.

## Porting to a new environment
To move to a new environment, you need to update the following files:
  + run/templates/<New_supercomputer_name>.sh
  + run/lib/manager.py
  + run/lib/supercomputers/<New_supercomputer_name>.py (newly added for new platform)
  + run/lib/supercomputers/supercomputers.json
  + run/jobs/<New_supercomputer_name>.json

### Add job template <New_supercomputer_name>.sh
Before creating the job template file, it is highly recommended to start with a naive shell script to find
the correct sets of compilers and environmental variables to run the executable on your environment.
The most importantly, we need device to device MPI communication on GPU environemnts. 

Once the working job script is prepared, the template job script should be created based on the sample job script.
```
cd run/templates
cp sub_Aquarius_A100_cuda.sh sub_<New_supercomputer_name>_<Device_name>_<Programming_Framework_name>.sh
```
The job template should depend on the scheduler, but in general it is necessary to make some variables (such as ```node```, ```ACCOUNT```, and so on) replacable. 
```bash
#!/bin/bash
#PJM -L "node={0[node]}"
#PJM -L "rscgrp={0[resource_group]}"
#PJM -L "elapse={0[time]}"
#PJM -s
#PJM -g {0[ACCOUNT]}
#PJM --mpi proc={0[mpi_procs]}

# Load modules
# default: module load gcc/8.3.1 cuda/11.2 ompi-cuda/4.1.1-11.2
module purge
{0[LOAD_MODULES]}

# environmental variables
# default: export UCX_MEMTYPE_CACHE=n; export UCX_IB_GPU_DIRECT_RDMA=no
{0[ENV_VARIABLES]}

echo start citylbm
mpiexec -machinefile $PJM_O_NODEINF -np $PJM_MPI_PROC -npernode {0[Processes_per_node]} ../bin/citylbm.A100_cuda \
{0[CITYLBM_OPTIONS]}> ../log/logfile.txt 2>&1

echo end citylbm
```

### Add <New_supercomputer_name>.py and update manager.py
You may copy an example as .
```
cd run/lib/supercomputers
cp _aquarius.py _<New_supercomputer_name>.py
```
Most of the script can be reused, but at least you need to modify the initialization part. 
```python3
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.machine_name = 'Aquarius'
        self.indent_params = 8
        self.default_modules = ['gcc/8.3.1',
                                'cuda/11.2',
                                'ompi-cuda/4.1.1-11.2',
                               ]
        self.default_envs = ['UCX_MEMTYPE_CACHE=n',
                             'UCX_IB_GPU_DIRECT_RDMA=no'] 
```
It is recommended to set default environmental variables and modules based on a workin job script.
Otherwise, these variables can be given from the input file.

Once the script dedicated to the <New_supercomputer_name> is made, 
it should be added to manager.py as follows
```python3
from .supercomputers._tsubame import Tsubame3
from .supercomputers._flow import Flow
from .supercomputers._sgi8600 import SGI8600
from .supercomputers._aquarius import Aquarius
from .supercomputers._<New_supercomputer_name> import <New_supercomputer_name>

def get_job_manager(name):
    SUPERCOMPUTERS = {
        'Tsubame3.0': Tsubame3,
        #'Fugaku': Fugaku,
        'SGI8600': SGI8600,
        'Flow': Flow,
        'Aquarius': Aquarius,
        '<New_supercomputer_name>': <New_supercomputer_name>,
        #'Odyssey': Odyssey,
    }
```

### Update Job class information
You need to add the Job class and device information of a new platform to ```run/lib/supercomputers/supercomputers.json.```
For example, the job class of `Aquarius` is given in the json file as
```json
{
"Aquarius": {
        "CPU_name": "Icelake",
        "CPUs_per_node": 2,
        "Cores_per_cpu": 36,
        "GPU_name": "A100",
        "GPUs_per_node": 8,
        "Job_classes": {
            "regular-a0": [
                1,
                2,
                "12:00:00"
            ],
            "regular-a1": [
                3,
                4,
                "12:00:00"
            ],
            "regular-a2": [
                5,
                8,
                "6:00:00"
            ]
        },
        "batch_command": "pjsub"
    },
}
```
The json file must include the CPU and GPU (if available) information. 
The Job class information in "Job_classes" has the structure of 
```
            "<job_class_name>": [
                <Minimum number of nodes in this job class>,
                <Maximum number of nodes in this job class>,
                "<Maximum wall time (in hh:mm:ss format) in this job class>",
            ],

```
.

### Preparing the sample input file <New_supercomputer_name>.json
The input file is a json file like
```json
{
    "JOB_settings": {
        "ACCOUNT": "gi37",
        "JOB_NAME": "oklahoma4m",
        "NB_PROCS": 8,
        "NB_THREADS": 1,
        "OUT_DIR": "/work/gr21/i18048/citylbm_base",
        "TIME": "10:00",
        "USE_GPUs": true
    },
    "citylbm_settings": {
        "cfr_flags": [
            1,
            1
        ],
        "cfr_steps": [
            750,
            750
        ],
        "domain_length": [
            4096,
            4096,
            2560
        ],
        "domain_min": [
            -2048,
            -2048,
            -8
        ],
        "n_ensemble_members": 1,
        "number_of_grid_point": [
            40,
            2
        ],
        "ofs_ensemble_idx": 0,
        "restart_flags_and_step": [
            0,
            99999
        ],
        "time_end": 24300.0,
        "velocity_lbm": [
            10.0,
            0.1
        ]
    }
}
```

Basically, the json file consists of two major categories ```JOB_settings``` and ```citylbm_settings```. 
See [the guide for input parameter settings](https://github.com/hasegawa-yuta-jaea/citylbm/blob/master/doc/Input-File-Guide.md) for detail.
1. ```JOB_settings```  
These are parameters for job script. 

1. ```citylbm_settings```  
These are physical and numerical parameters for citylbm.
