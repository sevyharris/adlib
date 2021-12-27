# should probable specify fcc111 for slab type

import sys
import yaml
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
from time import time


start = time()

settings_file = sys.argv[1]
with open(settings_file) as f:
    settings = yaml.load(f, Loader=yaml.FullLoader)

SLAB_DIR = settings['SLAB_DIR']
logfile = os.path.abspath(os.path.join(SLAB_DIR, 'ase.log'))

bulk_file = settings['BULK_FILE']
with open(bulk_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
initial_energy = traj[-1].get_potential_energy()
lattice_constant = traj[-1].cell[0][0]
print(f"Initial Energy: {initial_energy}")
print(f"Lattice constant: {lattice_constant}")


# Construct the slab + vacuum
vac = settings['vacuum_slab']
metal_slab = fcc111(
    settings['metal'],
    size=settings['slab_size'],
    vacuum=vac,
    a=lattice_constant
)

traj_file = os.path.abspath(os.path.join(SLAB_DIR, 'slab.traj'))
if os.path.exists(traj_file):
    traj = Trajectory(traj_file)
    if len(traj) > 0:
        adsorbate = traj[-1]

# Fix the bottom layer
bottom_layer = []
for i, pos in enumerate(metal_slab.get_positions()):
    if pos[-1] == vac:
        print(metal_slab[i])
        bottom_layer.append(metal_slab[i].index)

fix_bottom_layer = FixAtoms(indices=bottom_layer)
metal_slab.set_constraint(fix_bottom_layer)

espresso_settings = settings['slab_espresso_settings']

calc = Espresso(
    pseudopotentials=settings['pseudopotentials'],
    tstress=True,
    tprnfor=True,
    kpts=settings['kpts_slab'],
    pseudo_dir=settings['PSEUDOS_DIR'],
    input_data=espresso_settings,
    directory=SLAB_DIR
)

metal_slab.calc = calc
opt = BFGS(metal_slab, logfile=logfile, trajectory=traj_file)
opt.run(fmax=settings['forc_conv_thr_eVA'])

final_energy = metal_slab.get_potential_energy()
print(f"Final Energy: {final_energy}")

# copy the file
copyfile(os.path.abspath(os.path.join(SLAB_DIR, 'espresso.pwo')), settings['SLAB_FILE'])

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')
