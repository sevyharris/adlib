# should probable specify fcc111 for slab type

import sys
import os
from shutil import copyfile
from ase.build import bulk, fcc111
from ase.io import read, write
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms
from ase.io.espresso import read_espresso_out
from ase.io.trajectory import Trajectory
from ase.io.ulm import InvalidULMFileError
from time import time


start = time()

logfile = os.path.abspath(os.path.join(SLAB_DIR, 'ase.log'))

bulk_file = 'bulk.pwo'
with open(bulk_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
lattice_constant = traj[-1].cell[0][0]
print(f"Lattice constant: {lattice_constant}")

adsorbate_file = 'adsorbate.pwo'
with open(bulk_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
adsorbate = traj[-1]


# Construct the slab + vacuum
fmax = 0.001  # eV/A
vacuum = 15
height = 4.0
metal_slab = fcc111('Cu', size=(3, 3, 3), vacuum=vacuum, a=lattice_constant)

traj_file = 'slab.traj'
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            adsorbate = traj[-1]
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)

# Fix the bottom layer
bottom_layer = []
for i, pos in enumerate(metal_slab.get_positions()):
    if pos[-1] == vacuum:
        print(metal_slab[i])
        bottom_layer.append(metal_slab[i].index)

fix_bottom_layer = FixAtoms(indices=bottom_layer)
metal_slab.set_constraint(fix_bottom_layer)


# place the adsorbate
element_priority = ['C', 'O', 'H']  # where does N fit?
bond_atom_index = -1
for element in atom_priority:
    for atom in adsorbate:
        if atom.symbol == element:
            bond_atom = atom.index
            break
    if bond_atom_index > -1:
        break

# 'ontop', 'bridge', 'fcc', 'hcp'
# https://wiki.fysik.dtu.dk/ase/ase/build/surface.html#ase.build.add_adsorbate
add_adsorbate(metal_slab, adsorbate, height=height, position='ontop', mol_index=bond_atom_index)

write(f'initial_system.xyz', metal_slab)

espresso_settings = {
    'control': {
        'verbosity': 'high',
        'calculation': 'scf',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'degauss': 0.1,
        'ecutwfc': 100,
    },
}

pseudopotentials = {
    'H': 'H_ONCV_PBE-1.2.upf',
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'N': 'N_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'Pt': 'Pt_ONCV_PBE-1.2.upf',
}

pw_executable = os.environ['PW_EXECUTABLE']
calc = Espresso(
    command=f'{pw_executable} -in PREFIX.pwi > PREFIX.pwo'
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(4, 4, 1),
    pseudo_dir=os.environ['PSEUDO_DIR'],
    input_data=espresso_settings,
)

metal_slab.calc = calc
opt = BFGS(metal_slab, logfile=logfile, trajectory=traj_file)
opt.run(fmax=fmax)

final_energy = metal_slab.get_potential_energy()
print(f"Final Energy: {final_energy}")

# copy the file
copyfile('espresso.pwo', 'adsorbed.pwo')

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')
