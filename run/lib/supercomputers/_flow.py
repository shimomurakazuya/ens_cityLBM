from ._base import AbstractSupercomputer

class Flow(AbstractSupercomputer):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.machine_name = 'Flow'
        self.indent_params = 9
        
        # Compute how many nodes are required
        nb_nodes = 0
        nb_procs = self.job_dict['NB_PROCS']
        nb_threads = self.job_dict['NB_THREADS']
        
        try:
            self.json_data = self.json_data[self.machine_name]
        except:
            raise ValueError(f'Machine {self.machine_name} not registered to supercomputers.json')
        
        CPUs_per_node = self.json_data['CPUs_per_node']
        Cores_per_cpu = self.json_data['Cores_per_cpu']
        n_cores_per_node = CPUs_per_node * Cores_per_cpu
        nb_nodes = (nb_procs * nb_threads) // n_cores_per_node
        if (nb_procs * nb_threads) % n_cores_per_node != 0:
            raise ValueError('NB_PROCS * NB_THREADS in inputfile must be the multiple of CPU Cores_per_nodes: {n_cores_per_node}')
        template_file = 'templates/sub_Flow_A64FX_omp.sh'
        self.job_class, self.time_formatted = super()._get_job_class(nb_nodes, self.job_dict['TIME'])
        self.nb_nodes = nb_nodes
        
        ### Update job template
        parameters = {}
        parameters = self.job_dict.copy()
        parameters['node'] = self.nb_nodes
        parameters['time'] = self.time_formatted
        parameters['mpi_procs'] = nb_procs
        parameters['resource_group'] = self.job_class
        parameters['CITYLBM_OPTIONS'] = super()._options2str(self.citylbm_dict, default_indent=self.indent_params)
        
        with open(template_file, 'r') as f:
            template = f.read()
        template = template.format(parameters)
        
        job_script = self.exec_dir / 'job.sh'
        with open(job_script, 'w') as f:
            f.write(template)

        ### submission command
        args = []
        args.append('pjsub')
        args.append('job.sh')
        self.args = args
