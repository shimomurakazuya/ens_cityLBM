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
