import math
from ._base import AbstractSupercomputer

class FX700(AbstractSupercomputer):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.machine_name = 'FX700'
        self.default_indent = 8
        self.default_modules = ['openmpi4.1.1',
                               ]
        self.default_envs = ['OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}']

        try:
            self.json_data = self.json_data[self.machine_name]
        except:
            raise ValueError(f'Machine {self.machine_name} not registered to supercomputers.json')

        exec_command = self.json_data.get('mpi_command', 'mpirun')
        self.batch_command = self.json_data.get('batch_command', 'sbatch')

        # 
        self.default_dict = {
                             'modules': self.default_modules,
                             'envs': self.default_envs,
                             'exec_command': exec_command,
                             'app': 'citylbm.MI100_hip'
                            } # [TO DO] to be renamed to 'citylbm.exe' on every machine

        
        # Special treatment on FX700
        visible_devices = [1, 2, 3]
        visible_devices = [str(visible_device) for visible_device in visible_devices]
        visible_devices = 'ROCR_VISIBLE_DEVICES=' + ','.join(visible_devices)

        self.citylbm_settings = {}
        self.citylbm_settings['MPI_COMMAND'] = f'{visible_devices} {exec_command}'
        self.citylbm_settings['MPI_OPTIONS'] = ' -n ${SLURM_NTASKS}'
        self.citylbm_settings['REDIRECTION_COMMAND'] = '> ../log/logfile.txt 2>&1'

        if self.post_script:
            self.__init_post()
        else:
            self.__init_run()

    def __init_post(self):
        # Store the data in csv format
        raise NotImplementedError()

    def __init_run(self):
        # If something is wrong raise an Error
        self.__prepare(test_only=True)

        # Preparing directories
        self.__prepare()

        # submission command
        args = []
        args.append('sbatch')
        args.append('job.sh')
        self.args = args

    def __prepare(self, test_only=False):
        # Compute how many nodes are required
        nb_procs_list = self.job_dict['NB_PROCS']
 
        if type(nb_procs_list) is not list:
            # Then, production run
            nb_procs_list = [nb_procs_list]

        nb_procs_base = nb_procs_list[0]
        self.nb_nodes = []
        for nb_procs in nb_procs_list:
            # job_dict updated in the following function
            job_dict, template_file = super()._get_template_and_nodes(nb_procs=nb_procs)
            node = job_dict['node']
            self.nb_nodes.append(node)
 
            # Call this after _get_template_and_nodes()
            if not test_only:
                super()._prepare_citylbm_dirs(job_dict=job_dict)
 
            job_class, time_formatted = super()._get_job_class(node, job_dict['TIME'])
 
            ### Update citylbm settings for scalability measurement
            eps = 1.e-6
            scale = int(math.sqrt(nb_procs / nb_procs_base + eps))
            if nb_procs % nb_procs_base != 0:
                raise ValueError(f'NB_PROCS list in inputfile must be the multiple of 0th element: {nb_procs_base}')
 
            if not test_only:
                citylbm_dict = super()._update_citylbm_settings(citylbm_dict=self.citylbm_dict, scale=scale)
 
                ### Update job template
                ### Need to 
                parameters = job_dict.copy()
                parameters['resource_group'] = job_class
                parameters['time']           = time_formatted
                parameters['mpi_procs']      = nb_procs
                parameters['nodelist']       = 'amd1'

                ### Update job dict
                job_dict["EXECUTABLE_DIR"] = '../bin'
                job_dict['citylbm_settings'] = self.citylbm_settings.copy()

                parameters = super()._set_job_parameters(
                    parameters = parameters,
                    job_dict = job_dict,
                    default_job_dict = self.default_dict,
                    options_dict = citylbm_dict,
                    indent = self.default_indent
                )

                super()._generate_job_script(template_file, parameters)
