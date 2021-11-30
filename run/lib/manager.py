from .supercomputers._tsubame import Tsubame3
from .supercomputers._flow import Flow
from .supercomputers._sgi8600 import SGI8600
from .supercomputers._aquarius import Aquarius

def get_job_manager(name):
    SUPERCOMPUTERS = {
        'Tsubame3.0': Tsubame3,
        #'Fugaku': Fugaku,
        'SGI8600': SGI8600,
        'Flow': Flow,
        'Aquarius': Aquarius,
        #'Odyssey': Odyssey,
    }

    for n, supercomputer in SUPERCOMPUTERS.items():
        # Compare as lowercase
        if n.lower() == name.lower():
            return supercomputer

    raise ValueError(f'supercomputer {name} is not defined')
