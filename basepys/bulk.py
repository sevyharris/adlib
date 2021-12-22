# Module to compute the bulk lattice constant
# Sevy Harris
# 2021-12-03

import os
import sys
import yaml
from time import time
from ase.build import bulk
from ase.calculators.espresso import Espresso


start = time()
# read in the settings
settings_file = sys.argv[1]
with open(settings_file) as f:
    settings = yaml.load(f, Loader=yaml.FullLoader)

BULK_DIR = settings['BULK_DIR']
logfile = os.path.join(BULK_DIR, 'ase.log')

# N_CPUS = 4
# PW_EXEC = '/shared/centos7/q-e/q-e-qe-6.2.0-gcc7/bin/pw.x'
# input_file = os.path.join(BASEDIR, 's1_bulk', 'bulk.pwi')
# output_file = os.path.join(BASEDIR, 's1_bulk', 'bulk.pwo')
# os.environ['ASE_ESPRESSO_COMMAND'] = f'mpirun -np {N_CPUS} {PW_EXEC} -in {input_file} > {output_file}'
# os.environ['ASE_ESPRESSO_COMMAND'] = f'mpirun -np {N_CPUS} {PW_EXEC} -in PREFIX.pwi > PREFIX.pwo'

metal = settings['metal']
cu_bulk = bulk(settings['metal'],
               crystalstructure=settings['crystal_structure'],
               a=settings['lattice_constant_guess'],
               cubic=True
)


forc_conv_thr = settings['forc_conv_thr_eVA'] / 51.42208619083232
espresso_settings = {
    'control': {
        'calculation': 'vc-relax',
        'forc_conv_thr': forc_conv_thr,
    },
    'system': {
        'occupations': 'smearing',  # required for metals
        'degauss': 0.1,
        'ecutwfc': settings['ecutwfc'],  # 50
        'ecutrho': settings['ecutrho'],  # 500
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    },
    'cell': {
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'press_conv_thr': 0.5,
    }
}
if settings['dft_functional'].lower() != 'default':
    espresso_settings['system']['input_dft'] = settings['dft_functional']

bulk_calc = Espresso(pseudopotentials=settings['pseudopotentials'],
                     tstress=True,
                     tprnfor=True,
                     kpts=settings['kpts_bulk'],
                     pseudo_dir=settings['pseudos_dir'],
                     input_data=espresso_settings,
                     directory=BULK_DIR)
cu_bulk.calc = bulk_calc
cu_bulk.get_potential_energy()

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')
