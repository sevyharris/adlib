import os
import sys
import glob
from time import time
from ase.calculators.espresso import Espresso
from ase.optimize import BFGS
from ase.io.trajectory import Trajectory
from ase.io import read
from ase.io.espresso import read_espresso_out


start = time()
logfile = 'ase.log'


# Assume there is only one .xyz file in the folder
this_dir = os.path.dirname(os.path.abspath(__file__))
files = glob.glob(os.path.join(this_dir, '*.xyz'))
if len(files) > 1:
    print('multiple xyzs')


fmax = 0.01  # eV/A
vacuum = 10.0   # A
adsorbate_file = files[0]
adsorbate = read(adsorbate_file)
adsorbate.center(vacuum=vacuum)


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
        'ecutwfc': 500,
    },
}


pw_executable = os.environ['PW_EXECUTABLE']

pseudopotentials = {
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'N': 'N_ONCV_PBE-1.2.upf',
    'H': 'H_ONCV_PBE-1.2.upf',
}


# Restart if available
traj_file = 'ads.traj'
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            adsorbate = traj[-1]
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)

command = f'mpirun -np 16 {pw_executable} -in PREFIX.pwi > PREFIX.pwo'
print(command)

espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(1, 1, 1),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)


opt = BFGS(adsorbate, trajectory=traj_file, logfile=logfile)
adsorbate.calc = espresso
opt.run(fmax)


energy = 0
with open(os.path.join(this_dir, 'espresso.pwo'), 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    energy = traj[-1].get_potential_energy()


end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Energy: {energy} eV\n')
    f.write(f'Completed in {duration} seconds\n')

