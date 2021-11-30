import math
from ._base import AbstractSupercomputer

class SGI8600(AbstractSupercomputer):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.machine_name = 'SGI8600'
        self.indent_params = 7

        if not self.job_dict['USE_GPUs']:
            machine_name = 'SGI8600_CPU'
        else:
            machine_name = 'SGI8600'

        try:
            self.json_data = self.json_data[machine_name]
        except:
            raise ValueError(f'Machine {machine_name} not registered to supercomputers.json')

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

        ### Submission command
        args = []
        args.append('qsub')
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
                self._prepare_citylbm_dirs(job_dict=job_dict)

            job_class, time_formatted = super()._get_job_class(node, job_dict['TIME'])

            ### Update citylbm settings for scalability measurement
            eps = 1.e-6
            scale = int(math.sqrt(nb_procs / nb_procs_base + eps))
            if nb_procs % nb_procs_base != 0:
                raise ValueError(f'NB_PROCS list in inputfile must be the multiple of 0th element: {nb_procs_base}')
 
            if not test_only:
                citylbm_dict = super()._update_citylbm_settings(citylbm_dict=self.citylbm_dict, scale=scale)

                ### Update job template
                parameters = job_dict.copy()
                parameters['ACCOUNT']          = parameters['ACCOUNT'].upper() # in case account is given in lower case
                parameters['time']             = time_formatted
                parameters['mpi_procs']        = nb_procs
                parameters['resource_group']   = job_class
                parameters['ofs_ensemble_idx'] = citylbm_dict['ofs_ensemble_idx']
                parameters['CITYLBM_OPTIONS']  = super()._options2str(citylbm_dict, default_indent=self.indent_params)
                
                with open(template_file, 'r') as f:
                    template = f.read()
                template = template.format(parameters)
                
                job_script = self.exec_dir / 'job.sh'
                with open(job_script, 'w') as f:
                    f.write(template)
