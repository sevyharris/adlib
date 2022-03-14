import os
import sys
from time import time
from ase.calculators.espresso import Espresso
from ase.build import fcc111
from ase.constraints import FixAtoms
from ase.io.ulm import InvalidULMFileError
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from ase.io.espresso import read_espresso_out


start = time()
logfile = 'ase.log'
this_dir = os.path.dirname(os.path.abspath(__file__))


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
        'ecutwfc': 100,
    },
}


fmax = 0.01
vacuum = 10.0
slab = fcc111('Cu', size=(3, 3, 4), vacuum=vacuum, a=3.6)


# Fix the bottom two layers
bottom_layer = []
second_layer = []
fixed_indices = []
z_values = list(set([pos[2] for pos in slab.get_positions()]))
z_values.sort()


for i, pos in enumerate(slab.get_positions()):
    if pos[2] == z_values[0]:
       bottom_layer.append(slab[i].index)
    if pos[2] == z_values[1]:
       second_layer.append(slab[i].index)
fixed_indicies = bottom_layer + second_layer
fix_bottom_layers = FixAtoms(indices=fixed_indicies)
slab.set_constraint(fix_bottom_layers)


# Restart if available
traj_file = 'slab.traj'
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            slab = traj[-1]
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)


pw_executable = os.environ['PW_EXECUTABLE']

pseudopotentials = {
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'N': 'N_ONCV_PBE-1.2.upf',
    'H': 'H_ONCV_PBE-1.2.upf',
    'Pt': 'Pt_ONCV_PBE-1.2.upf',
    'Pd': 'Pd_ONCV_PBE-1.2.upf',
}

command = f'mpirun -np 32 {pw_executable} -in PREFIX.pwi > PREFIX.pwo'
print(command)

espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(7, 7, 7),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

slab.calc = espresso
opt = BFGS(slab, logfile=logfile, trajectory='slab.traj')
opt.run(fmax=fmax)

# Read the energy back in
energy = 0
with open(os.path.join(this_dir, 'espresso.pwo'), 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    energy = traj[-1].get_potential_energy()


end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Energy: {energy} eV\n')
    f.write(f'Completed in {duration} seconds\n')

