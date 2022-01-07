# Module to compute the bulk lattice constant
# Sevy Harris
# 2021-12-03

import os
import shutil
import sys
import yaml
from shutil import copyfile
from time import time
from ase.calculators.espresso import Espresso
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.io.trajectory import Trajectory
from ase.io.ulm import InvalidULMFileError


start = time()
# read in the settings
settings_file = sys.argv[1]
with open(settings_file) as f:
    settings = yaml.load(f, Loader=yaml.FullLoader)

ADS_DIR = settings['ADS_DIR']
logfile = os.path.abspath(os.path.join(ADS_DIR, 'ase.log'))


espresso_settings = settings['ads_espresso_settings']
# remove everything that has 'default' as the value
for category_key in espresso_settings.keys():
    keys_to_remove = []
    for setting_key in espresso_settings[category_key].keys():
        if espresso_settings[category_key][setting_key] == 'default':
            keys_to_remove.append(setting_key)
    for key in keys_to_remove:
        espresso_settings[category_key].pop(key)


adsorbate = read(settings['adsorbate_file'])
adsorbate.center(vacuum=settings['vacuum_ads'])

# Restart if available
traj_file = os.path.abspath(os.path.join(ADS_DIR, 'ads.traj'))
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            adsorbate = traj[-1]
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)

ads_calc = Espresso(
    pseudopotentials=settings['pseudopotentials'],
    tstress=True,
    tprnfor=True,
    kpts=settings['kpts_ads'],
    pseudo_dir=settings['PSEUDOS_DIR'],
    input_data=espresso_settings,
    directory=ADS_DIR
)
adsorbate.calc = ads_calc

opt = BFGS(adsorbate, trajectory=traj_file, logfile=logfile)
opt.run(fmax=settings['forc_conv_thr_eVA'])

# copy the file
copyfile(os.path.abspath(os.path.join(ADS_DIR, 'espresso.pwo')), settings['BULK_FILE'])

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')
