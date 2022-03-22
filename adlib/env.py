import os
import warnings


# TODO add max_cpu settings for personal computer vs. toastie
discovery_env = {
    'name': 'DISCOVERY',
    'max_cpus': 64,
}
single_node_env = {
    'name': 'SINGLE_NODE',
    'max_cpus': 4,
}

# NOTE all environments require openmpi
valid_environments = [
    'SINGLE_NODE',  # personal computer -- default
    'DISCOVERY',    # SLURM cluster
    'THETA',        # COBALT cluster
]
environments = {
    'SINGLE_NODE': single_node_env,
    'DISCOVERY': discovery_env,
}


def load_environment():
    try:
        environment = os.environ['ADLIB_ENV']
        if environment.upper() not in valid_environments:
            raise NotImplementedError(f'Environment {environment.upper()} not in valid environments: {valid_environments}')
        return environment.upper()
    except KeyError:
        warnings.warn(f'No ADLIB_ENV speficied. Using SINGLE_NODE...')
        # raise('No ADLIB_ENV specified. Set your environment with "export ADLIB_ENV=SINGLE_NODE"')
        return valid_environments[0]
