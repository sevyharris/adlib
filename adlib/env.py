import os
import warnings


# NOTE all environments require openmpi
valid_environments = [
    'SINGLE_NODE',  # personal computer -- default
    'DISCOVERY',    # SLURM cluster
    'THETA',        # COBALT cluster
]


def load_environment():
    try:
        environment = os.environ['ADLIB_ENV']
        if environment.upper not in valid_environments:
            raise NotImplementedError(f'Environment {environment.upper} not in valid environments: {valid_environments}')
        return environment.upper
    except KeyError:
        warnings.warn(f'No ADLIB_ENV speficied. Using SINGLE_NODE...')
        # raise('No ADLIB_ENV specified. Set your environment with "export ADLIB_ENV=SINGLE_NODE"')
        return valid_environments[0]
