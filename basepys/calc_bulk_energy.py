import os
import sys
from time import time
from ase.calculators.espresso import Espresso
from ase.build import bulk


start = time()
logfile = 'ase.log'


espresso_settings = {
    'control': {
        'verbosity': 'high',
        'calculation': 'scf',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'smearing': 'mv',
        'degauss': 0.1,
        'ecutwfc': 1000,
    },
}


cu_bulk = bulk('Cu', crystalstructure='fcc', a=3.5500000000000003, cubic=True)

pw_executable = os.environ['PW_EXECUTABLE']

pseudopotentials = {
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'N': 'N_ONCV_PBE-1.2.upf',
    'H': 'H_ONCV_PBE-1.2.upf',
}

command = f'mpirun -np 16 {pw_executable} -in PREFIX.pwi > PREFIX.pwo'
print(command)

espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(9, 9, 9),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

cu_bulk.calc = espresso
energy = cu_bulk.get_potential_energy()


end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Energy: {energy} eV\n')
    f.write(f'Completed in {duration} seconds\n')

