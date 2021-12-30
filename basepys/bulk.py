# Module to compute the bulk lattice constant
# Sevy Harris
# 2021-12-03

import os
from shutil import copyfile
import sys
import yaml
from time import time
from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
import socket

start = time()
# read in the settings
settings_file = sys.argv[1]
with open(settings_file) as f:
    settings = yaml.load(f, Loader=yaml.FullLoader)

BULK_DIR = settings['BULK_DIR']
logfile = os.path.abspath(os.path.join(BULK_DIR, 'ase.log'))

# N_CPUS = 4
# PW_EXEC = '/shared/centos7/q-e/q-e-qe-6.2.0-gcc7/bin/pw.x'
# input_file = os.path.join(BASEDIR, 's1_bulk', 'bulk.pwi')
# output_file = os.path.join(BASEDIR, 's1_bulk', 'bulk.pwo')
# os.environ['ASE_ESPRESSO_COMMAND'] = f'mpirun -np {N_CPUS} {PW_EXEC} -in {input_file} > {output_file}'
# os.environ['ASE_ESPRESSO_COMMAND'] = f'mpirun -np {N_CPUS} {PW_EXEC} -in PREFIX.pwi > PREFIX.pwo'

metal = settings['metal']
cu_bulk = bulk(
    settings['metal'],
    crystalstructure=settings['crystal_structure'],
    a=settings['lattice_constant_guess'],
    cubic=True
)

espresso_settings = settings['bulk_espresso_settings']
espresso_settings['control']['forc_conv_thr'] = settings['forc_conv_thr_eVA'] / 51.42208619083232

# remove everything that has 'default' as the value
for category_key in espresso_settings.keys():
    keys_to_remove = []
    for setting_key in espresso_settings[category_key].keys():
        if espresso_settings[category_key][setting_key] == 'default':
            keys_to_remove.append(setting_key)
    for key in keys_to_remove:
        espresso_settings[category_key].pop(key)


pw_executable = settings['pw_executable']
port = 30141
hostname = socket.gethostname()

espresso = Espresso(
    command=f'{pw_executable} -in PREFIX.pwi --ipi {hostname}:{port} > PREFIX.pwo',
    pseudopotentials=settings['pseudopotentials'],
    tstress=True,
    tprnfor=True,
    kpts=settings['kpts_bulk'],
    pseudo_dir=settings['PSEUDOS_DIR'],
    input_data=espresso_settings,
    directory=BULK_DIR
)

with SocketIOCalculator(espresso, log=sys.stdout) as calc:
    cu_bulk.calc = calc
    cu_bulk.get_potential_energy()


# copy the file
copyfile(os.path.abspath(os.path.join(BULK_DIR, 'espresso.pwo')), settings['BULK_FILE'])

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')
