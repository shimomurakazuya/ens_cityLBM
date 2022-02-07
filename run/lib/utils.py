import subprocess
from distutils.util import strtobool

def str2num(datum):
    """ Convert string to integer or float"""

    try:
        return int(datum)
    except:
        try:
            return float(datum)
        except:
            try:
                return strtobool(datum)
            except:
                return datum

def list2str(var_list, prefix, suffix='\n'):
    """
    ex. var_list = ['gcc/8.3.1', 'cuda/11.2', 'ompi-cuda/4.1.1-11.2']
        prefix = 'load module'
        return 'load module gcc/8.3.1\nload module cuda/11.2\nload module ompi-cuda/4.1.1-11.2\n'
    """

    if type(var_list) is not list:
        var_list = [var_list]

    if type(prefix) is list:
        prefix = ' '.join(prefix)

    tmp_str = ''
    for value in var_list:
        tmp_str += f'{prefix} {value}{suffix}'

    return tmp_str

def dict2str(dict, indent, separator='--', add_new_line_in_the_end=False):
    """
    Convert dict to a single string containing arguments of the app
    """

    tmp_str = ''
    max_len = len(max(dict.keys(), key=lambda x: len(x))) + 5 # add some space for indenting

    last_element = list(dict.keys())[-1]
    for option, values in dict.items():
        if type(values) is not list:
            values = [values]

        # Corresponding to the argument name, list "-n_ensemble_members"
        tmp_str += f'{"":<{indent}}{separator}{option: <{max_len}}'

        # Corresponding to values, like 100 or -2048 -2048 -8
        for value in values:
            tmp_str += f'\t {value}'

        # New line
        if add_new_line_in_the_end or (option is not last_element):
            tmp_str += ' \\\n'
    return tmp_str

def get_meta():
    """
    Get git version, date, user, and hostname
    """
    meta = ''

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

def to_time_format(time_str, separator=':'):
    """
    Parameters
    ----------
    time_str : string
        Wall time in string given in the following form
        1. 'hh:mm:ss', 2. 'mm:ss', 3. 'mm'

    Returns
    -------
    time_str : string
        Wall time in string given in the following form
        'hh:mm:ss'
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

def less_eq_time(time_a, time_b, separator=':'):
    """
    Compare the time stored in 'hh:mm:ss' format
    """

    def in_seconds(time_str):
        time_list = time_str.split(separator)
        hh, mm, ss = time_list
        seconds = str2num(hh) * 3600 + str2num(mm) * 60 + str2num(ss)
        return seconds

    return in_seconds(time_a) <= in_seconds(time_b)
