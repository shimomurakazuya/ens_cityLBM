"""
Usage

Step0. (Porting to a new environment) Update the json script in lib, create a template job script file in template and modify lib/manager.py
Step1. Set environmental variable export SUPERCOMPUTER=<SGI8600 (default), Tsubame3.0, FlowType1, FlowType2, Fugaku, Aquarius, Odyssey>
Step2. Locate input file <your_input_file> at <your_input_file_dir> (default=./jobs)
Step3. Run python script as
    python sub.py <your_input_file> --dir <your_input_file_dir> (default=./jobs)
Step4. Then, the symbolic link to the execution directory will be created
"""

__author__     = 'Yuuichi ASAHI'
__date__       = '2022/01/27'
__version__    = '1.2'
__maintainer__ = 'Yuuichi ASAHI'
__email__      = 'asahi.yuichi@jaea.go.jp'
__status__     = 'Production'

import os
from lib.parser import parse
from lib.manager import get_job_manager

if __name__ == '__main__':
    args = parse()
    inputfile = args.inputfile
    inputfile_dir = args.dir
    verbose = args.verbose
    create_symlink = args.create_symlink
    machine_name = os.getenv('SUPERCOMPUTER', 'SGI8600')

    job_manager = get_job_manager(name=machine_name)(inputfile=inputfile, dirname=inputfile_dir, create_symlink=create_symlink)
    job_manager.submit(verbose=verbose)
