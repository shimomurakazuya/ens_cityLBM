import math
from ._base import AbstractSupercomputer

class Aquarius(AbstractSupercomputer):
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
        
        try:
            self.json_data = self.json_data[self.machine_name]
        except:
            raise ValueError(f'Machine {self.machine_name} not registered to supercomputers.json')

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

        ### submission command
        args = []
        args.append('pjsub')
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

            # Special treatment on aquarius and odyssey
            # remove the last letter since actual job class is regular-a rather than regular-a0
            job_class = job_class[:-1]
 
            ### Update citylbm settings for scalability measurement
            eps = 1.e-6
            scale = int(math.sqrt(nb_procs / nb_procs_base + eps))
            if nb_procs % nb_procs_base != 0:
                raise ValueError(f'NB_PROCS list in inputfile must be the multiple of 0th element: {nb_procs_base}')
 
            if not test_only:
                citylbm_dict = super()._update_citylbm_settings(citylbm_dict=self.citylbm_dict, scale=scale)
 
                ### Update job template
                parameters = job_dict.copy()
                parameters['resource_group']  = job_class
                parameters['time']            = time_formatted
                parameters['mpi_procs']       = nb_procs
                parameters['CITYLBM_OPTIONS'] = super()._options2str(citylbm_dict, default_indent=self.indent_params)

                # Adding modules and variables
                def get_values(key, default_values):
                    if key in self.job_dict:
                        return self.job_dict[key]
                    else:
                        return default_values

                modules = get_values(key='modules', default_values=self.default_modules)
                envs    = get_values(key='envs',    default_values=self.default_envs)

                parameters['LOAD_MODULES']    = super()._list2str(modules, prefix='module load')
                parameters['ENV_VARIABLES']   = super()._list2str(envs,    prefix='export')
 
                with open(template_file, 'r') as f:
                    template = f.read()
                template = template.format(parameters)
 
                job_script = self.exec_dir / 'job.sh'
                with open(job_script, 'w') as f:
                    f.write(template)
