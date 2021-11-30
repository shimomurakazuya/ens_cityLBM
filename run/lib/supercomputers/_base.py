"""Supercomputer base class
"""

import json
import os
import subprocess
import pathlib
import shutil
from ..parser import parse_from_file
from ..utils import str2num

class AbstractSupercomputer:
    def __init__(self, **kwargs):
        self.machine_name = None
        self.inputfile    = None
        self.exec_dirs    = [] # This empty list is updated and used for job submission stage
        
        allowed_kwargs = {
                          'inputfile',
                          'dirname',
                          'json_file',
                          'create_symlink',
                          'post_script',
                         }
        
        for kwarg in kwargs:
            if kwarg not in allowed_kwargs:
                raise TypeError('Keyword argument not understood: ', kwarg)
           
        inputfile = kwargs.get('inputfile')
        if not inputfile:
            IOError(f'inputfile {inputfile} is not specified')
        dirname = kwargs.get('dirname', 'jobs')
        self.post_script = kwargs.get('post_script', False)
        self.create_symlink = kwargs.get('create_symlink', True)
        self.json_file = kwargs.get('json_file', 'lib/supercomputers/supercomputers.json')
        var_dict = parse_from_file(dirname=dirname, filename=inputfile)
        self.job_dict = var_dict['JOB_settings']
        self.citylbm_dict = var_dict['citylbm_settings']
        self.dirname = dirname
        self.inputfile = inputfile
        
        with open(self.json_file, 'r') as f:
            self.json_data = json.load(f)

    def submit(self, verbose=True):
        for node, exec_dir in zip(self.nb_nodes, self.exec_dirs):
            os.chdir(exec_dir)
            if verbose:
                print(f'submitting {node} nodes job at {exec_dir}')
            subprocess.run(args=self.args)

    def __to_time_format(self, time_str, separator=':'):
        """
        Input:
            1. hh:mm:ss
            2. mm:ss
            3. mm
        Output:
            hh:mm:ss
        """
        time_list = time_str.split(separator)
 
        hh, mm, ss = "00", "00", "00"
        if len(time_list) == 3:
            hh, mm, ss = time_list
        elif len(time_list) == 2:
            mm, ss = time_list
        elif len(time_list) == 1:
            mm = time_list
        else:
            raise ValueError('TIME in the inputfile must be given in either hh:mm:ss, mm:ss, or mm format')
 
        return f'{hh}:{mm}:{ss}'

    def __less_eq_time(self, time_a, time_b, separator=':'):
        """
        Compare the time stored in hh:mm:ss format
        """
             
        def in_seconds(time_str):
            time_list = time_str.split(separator)
            hh, mm, ss = time_list
            seconds = str2num(hh) * 3600 + str2num(mm) * 60 + str2num(ss)
            return seconds
            
        return in_seconds(time_a) <= in_seconds(time_b)

    def  __get_meta(self):
        """
        Get git version, date, user, and hostname
        """
        meta = ""
 
        args_dict = {
                     'Git version       : ': ['git', 'rev-parse', 'master'],
                     'Date              : ': ['date'],
                     'User              : ': ['whoami'],
                     'Host name         : ': ['hostname'],
                     'Current directory : ': ['pwd'],
                    }
 
        for meta_name, args in args_dict.items():
            stdout = subprocess.run(args=args, encoding='utf-8', stdout=subprocess.PIPE).stdout
            meta += (meta_name + stdout)
 
        meta = meta[:-1] # removing the final '\n'
        return meta

    def _prepare_citylbm_dirs(self, job_dict):
        base_dir = pathlib.Path( job_dict['OUT_DIR'] )
        nb_procs = job_dict['NB_PROCS']
        if not base_dir.exists():
            base_dir.mkdir(parents=True)
        
        # Parent dir: base_dir/<JOB_NAME>_<DEVICE_NAME>_<NB_NODES>_<MPI_procs>_<SCALING>
        wk_dir = f'{job_dict["JOB_NAME"]}_{job_dict["DEVICE"]}_Nodes{job_dict["node"]}_MPI{nb_procs}'
        if job_dict['SCALING'] != 'None':
            wk_dir += f'_{job_dict["SCALING"]}'

        citylbm_dir = base_dir / wk_dir
        if not citylbm_dir.exists():
            citylbm_dir.mkdir(parents=True)
        print(f'--- Start preparing directories at {citylbm_dir} ---')
        
        # Child dirs
        dirs = ['run', 'bin', 'io', 'log']
        for dir in dirs:
            child_dir = citylbm_dir / dir
            if not child_dir.exists():
                child_dir.mkdir(parents=True)
        
        self.exec_dir = citylbm_dir / 'run'

        # Keep exec_dirs for job submission
        self.exec_dirs.append(self.exec_dir)

        # Copy input json file
        input_json = f'{self.dirname}/{self.inputfile}'
        shutil.copy(input_json, self.exec_dir)
        
        # Copy inputfiles
        io_dir = pathlib.Path('../io')
        for io_subdir in io_dir.iterdir():
            subdir = str(io_subdir).replace('../io/', '')
            if not (citylbm_dir / 'io' / subdir).exists():
                shutil.copytree(io_subdir, citylbm_dir / 'io' / subdir)
            else:
                print(f'inputdir {subdir} already exists')
        
        # Copy executable
        bin_dir = pathlib.Path('../bin')
        executable = list(bin_dir.glob('citylbm.*'))
        if not executable:
            raise FileNotFoundError('citylbm executable does not exist at bin')
        
        for ex in executable:
            shutil.copy(ex, citylbm_dir / 'bin')
        
        print(f'--- directories are ready at {citylbm_dir} ---')
        symdir = wk_dir
        if self.create_symlink and not os.path.exists(symdir):
            os.symlink(citylbm_dir, symdir)
            print(f'--- symbolic_link to {citylbm_dir}: {symdir} ---')

        # Save Meta data
        meta = self.__get_meta()
        meta_file = self.exec_dir / 'meta.txt'
        with open(meta_file, 'w') as f:
            f.write(meta)

    # Methods called from subclasses
    def _get_job_class(self, nb_nodes, time_str):
        env_dict = self.json_data
        job_classes = env_dict['Job_classes']
        
        job_class = None
        time_formatted = None

        # First format the time
        time_str = self.__to_time_format(time_str)
        
        # Find a job class that satisfies the node requirements
        max_nodes_in_job_class = 100000 # some large number
        for key, value in job_classes.items():
            min_nodes, max_nodes, max_time = value
        
            if min_nodes <= nb_nodes <= max_nodes and max_nodes < max_nodes_in_job_class:
                max_nodes_in_job_class = max_nodes
                # Use this job class
                if self.__less_eq_time(time_str, max_time):
                    job_class, time_formatted = key, time_str
                else:
                    print(f'{key} may be avilable if wall_time is reduced from {time_str} to {max_time}')
                    #break
                    #job_class, time_formatted = key, max_time
                        
        if not job_class:
            raise ValueError(f'No queue batch available for this node number {nb_nodes}')
                               
        return job_class, time_formatted
                                 
    def _options2str(self, dict, default_indent):
        """
        Convert options in dict to a single string containg options for citylbm
        """
        
        tmp_str = ''
        max_len = len(max(dict.keys(), key=lambda x: len(x))) + 5 # add some space for indenting
        for option, values in dict.items():
            if type(values) is not list:
                values = [values]
            # Corresponding to the argument name, like "-n_ensemble_members"
            tmp_str += f'{"":<{default_indent}}-{option: <{max_len}}'
            
            # Corresponding to the value, like 100 or -2048  -2048  -8
            for value in values:
                tmp_str += f'\t {value}'
                  
            # New line
            tmp_str += ' \\\n'
        return tmp_str

    def _list2str(self, var_list, prefix):
        """
        ex. var_list = ['gcc/8.3.1', 'cuda/11.2', 'ompi-cuda/4.1.1-11.2']
            prefix = 'load module'

            return 'load module gcc/8.3.1\nload module cuda/11.2\nload module ompi-cuda/4.1.1-11.2\n'
        """

        tmp_str = ''
        for value in var_list:
            tmp_str += f'{prefix} {value}\n'

        return tmp_str

    def _get_template_and_nodes(self, nb_procs):
        job_dict = self.job_dict.copy()
        nb_threads = job_dict['NB_THREADS']
        nb_cores = nb_threads * nb_procs

        CPUs_per_node = self.json_data['CPUs_per_node']
        Cores_per_cpu = self.json_data['Cores_per_cpu']
        nb_cores_per_node = CPUs_per_node * Cores_per_cpu

        if job_dict['USE_GPUs']:
            # GPU job
            try:
                device = self.json_data['GPU_name']
            except:
                raise ValueError(f'GPU information missing for {self.machine_name} in supercomputers.json')

            GPUs_per_node = self.json_data['GPUs_per_node']

            if nb_procs <= GPUs_per_node:
                # Single node GPU job
                nb_nodes = 1
                job_dict['Processes_per_node'] = nb_procs
            else:
                nb_nodes = nb_procs // GPUs_per_node
                if nb_procs % GPUs_per_node != 0:
                    raise ValueError(f'NB_PROCS in inputfile must be the multiple of GPUs_per_nodes: {GPUs_per_node}')
                job_dict['Processes_per_node'] = GPUs_per_node

            nb_cores = job_dict['Processes_per_node'] * nb_threads
            if nb_cores > nb_cores_per_node:
                raise ValueError(f'GPUs_per_node * NB_THREADS in inputfile must be smaller than CPU Cores_per_nodes: {nb_cores_per_node}')

            self.citylbm_dict['gpu_per_node'] = job_dict['Processes_per_node']
            job_dict['nb_cores'] = nb_cores
            template_file = f'templates/sub_{self.machine_name}_{device}_cuda.sh'
        else:
            # CPU job
            try:
                device = self.json_data['CPU_name']
            except:
                raise ValueError(f'CPU information missing for {self.machine_name} in supercomputers.json')

            if nb_cores <= nb_cores_per_node:
                # Single node CPU job
                nb_nodes = 1
                job_dict['Processes_per_node'] = nb_procs
            else:
                nb_nodes = nb_cores // nb_cores_per_node
                if nb_cores % nb_cores_per_node != 0:
                    raise ValueError(f'NB_PROCS * NB_THREADS in inputfile must be the multiple of CPU Cores_per_nodes: {nb_cores_per_node}')
                job_dict['Processes_per_node'] = nb_procs // nb_nodes

            template_file = f'templates/sub_{self.machine_name}_{device}_omp.sh'

        job_dict['DEVICE'] = device
        job_dict['node']   = nb_nodes
        job_dict['NB_PROCS'] = nb_procs

        return job_dict, template_file

    def _update_citylbm_settings(self, citylbm_dict, scale):
        """
        Change problem size for weak scaling measurement
        """
 
        tmp_citylbm_dict = citylbm_dict.copy()
        if self.job_dict['SCALING'] == 'Weak':
            # just enlarge the domain size
            lx, ly, lz          = citylbm_dict['domain_length']
            x_min, y_min, z_min = citylbm_dict['domain_min']
 
            scaleup = lambda x: int(x * scale)
            lx, ly       = scaleup(lx), scaleup(ly)
            x_min, y_min = scaleup(x_min), scaleup(y_min)
 
            tmp_citylbm_dict['domain_length'] = [lx, ly, lz]
            tmp_citylbm_dict['domain_min']    = [x_min, y_min, z_min]
 
        return tmp_citylbm_dict
