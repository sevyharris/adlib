# Module to compute the bulk lattice constant
# Sevy Harris
# 2021-12-03

import os
from shutil import copyfile
import sys
from time import time
from ase.build import bulk
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
import socket

start = time()

logfile = 'ase.log'


metal = settings['metal']
cu_bulk = bulk('Pt', 'fcc', a=4.0, cubic=True)
espresso_settings = {
    'control': {
        'verbosity': 'high',
        # 'disk_io': 'none',
        'calculation': 'vc-relax',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'degauss': 0.1,
        'ecutwfc': 500,
        # 'ecutrho': 500,
    },
    'ions': {
        'ion_dynamics': 'bfgs',
    },
    'cell': {
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'press_conv_thr': 0.005,
    }
}

pseudopotentials = {
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'Pt': 'Pt_ONCV_PBE-1.2.upf',
}

pw_executable = os.environ['PW_EXECUTABLE']
espresso = Espresso(
    command=f'mpirun -np 64 {pw_executable} -nimage 4 -npool 4 -in PREFIX.pwi > PREFIX.pwo',
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(16, 16, 16),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
    directory=BULK_DIR
)

cu_bulk.calc = calc
cu_bulk.get_potential_energy()


# copy the file so we don't accidentally overwrite it
copyfile('espresso.pwo', 'bulk.pwo')

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')
