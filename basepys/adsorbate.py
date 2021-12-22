# Module to compute the bulk lattice constant
# Sevy Harris
# 2021-12-03

import os
import sys
import yaml
from time import time
from ase.calculators.espresso import Espresso
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS


start = time()
# read in the settings
settings_file = sys.argv[1]
with open(settings_file) as f:
    settings = yaml.load(f, Loader=yaml.FullLoader)

ADS_DIR = settings['ADS_DIR']
logfile = os.path.join(ADS_DIR, 'ase.log')

espresso_settings = {
    'control': {
        'calculation': 'scf',
    },
    'system': {
        'ecutwfc': 50,
        'ecutrho': 500
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    },
}

if settings['dft_functional'].lower() != 'default':
    espresso_settings['system']['input_dft'] = settings['dft_functional']

adsorbate = read(settings['adsorbate_file'])
adsorbate.center(vacuum=settings['vacuum_ads'])

ads_calc = Espresso(pseudopotentials=settings['pseudopotentials'],
                    tstress=True,
                    tprnfor=True,
                    kpts=settings['kpts_ads'],
                    pseudo_dir=settings['pseudos_dir'],
                    input_data=espresso_settings,
                    directory=ADS_DIR)
adsorbate.calc = ads_calc

opt = BFGS(adsorbate, trajectory='ads.traj', logfile=logfile)
opt.run(fmax=settings['forc_conv_thr_eVA'])

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')
