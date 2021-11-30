import argparse
import json

def parse():
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', \
                        action='store', \
                        nargs=None, \
                        const=None, \
                        default=None, \
                        type=str, \
                        choices=None, \
                        help='CityLBM input file', \
                        metavar=None)
    parser.add_argument('--dir', \
                        action='store', \
                        nargs='?', \
                        const=None, \
                        default='./jobs', \
                        type=str, \
                        choices=None, \
                        help='Path for an input file (default: ./jobs)', \
                        metavar=None)
    parser.add_argument('--verbose', \
                        action='store_false', \
                        default=True, \
                        help='verbose job info or not')
    parser.add_argument('--create_symlink', \
                        action='store_false', \
                        default=True, \
                        help='Create symbolic link or not')

    return parser.parse_args()

def sanity_check(json_data, filename):
    """
    Check the inputfile includes setting parameters
    """
    try:
        job_dict = json_data['JOB_settings']
        citylbm_dict = json_data['citylbm_settings']
    except:
        raise IOError(f'JOB_settings or citylbm_settings is missing in the json input file {filename}')
                       
    required_job_settings = ['ACCOUNT', 'NB_PROCS', 'NB_THREADS', 'OUT_DIR', 'USE_GPUs']
    required_citylbm_settings = ['domain_length', 'domain_min', 'number_of_grid_point', 'time_end']
                          
    def check_missing_settings(required_list, dict):
        if not all(elem in dict.keys() for elem in required_list):
            raise IOError(f'one of {required_list} is missing in the json input file {filename}')
                                      
    check_missing_settings(required_job_settings, job_dict)
    check_missing_settings(required_citylbm_settings, citylbm_dict)

    # Set default parameters
    # If values are given, given values are used
    def set_default_values(dict, key, default_value):
        dict[key] = dict.get(key, default_value)
        return dict

    job_default_dict = {'JOB_NAME': 'oklahoma', 'TIME': '10:00:00', 'SCALING': 'None'}
    for key, value in job_default_dict.items():
        job_dict = set_default_values(job_dict, key, value)

    citylbm_default_dict = {
                            'gpu_per_node': 99999, # some large number, set as correct number if USE_GPUs is True
                            'cfr_flags': [1, 1],
                            'cfr_steps': [50, 50],
                            'n_ensemble_members': 1,
                            'ofs_ensemble_idx': 0,
                            'restart_flags_and_step': [0, 99999],
                            'velocity_lbm': [10.0, 0.1],
                           }

    for key, value in citylbm_default_dict.items():
        citylbm_dict = set_default_values(citylbm_dict, key, value)

    # Update original dict
    json_data['JOB_settings']     = job_dict
    json_data['citylbm_settings'] = citylbm_dict

    return json_data

def parse_from_file(dirname, filename):
    path = f'{dirname}/{filename}'
    with open(path, 'r') as f:
        json_data = json.load(f)

    json_data = sanity_check(json_data, path)
    return json_data
