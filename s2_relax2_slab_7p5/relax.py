# [1] https://doi.org/10.1063/1.4994149

from ase.build import bulk, fcc111
from ase.io import read, write
from ase import Atoms
from ase.optimize import LBFGS
from ase.calculators.espresso import Espresso
from ase.constraints import FixAtoms
from ase.io.espresso import read_espresso_out
from time import time



start = time()
logfile = 'ase.log'

# read in the results from the bulk cell relaxation
prev_run = 'run1.pwo'
with open(prev_run, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    cu_slab = traj[-1]

# bulk_file = '/home/harris.se/dft_adsorption/s1_relax_bulk/espresso.pwo'
#
# with open(bulk_file, 'r') as f:
#     traj = list(read_espresso_out(f, index=slice(None)))

# initial_energy = traj[-1].get_potential_energy()
# lattice_constant = traj[-1].cell[0][0]

# print(f"Initial Energy: {initial_energy}")
# print(f"Lattice constant: {lattice_constant}")


# Construct the slab + vacuum
vac = 7.5
#cu_slab = fcc111('Cu', size=(3, 3, 3), vacuum=vac, a=lattice_constant)

# Fix the bottom layer
bottom_layer = []
for i, pos in enumerate(cu_slab.get_positions()):
    if pos[-1] == vac:
        print(cu_slab[i])
        bottom_layer.append(cu_slab[i].index)

fix_bottom_layer = FixAtoms(indices=bottom_layer)
cu_slab.set_constraint(fix_bottom_layer)
write('slab_init.xyz', cu_slab)


pseudopotentials = {'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF'}
input_settings = {
    'control': {
        'calculation': 'relax',
        'forc_conv_thr': 0.001
    },
    'system': {
        'occupations': 'smearing',  # required for metals
        'degauss': 0.1,
        'ecutwfc': 50,
        'ecutrho': 500
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    },
}

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(4, 4, 1),
                pseudo_dir='/home/harris.se/espresso/pseudos/',
                input_data=input_settings)

cu_slab.calc = calc
opt = LBFGS(cu_slab, logfile=logfile, trajectory='qn.traj')
opt.run(fmax=0.001)


final_energy = cu_slab.get_potential_energy()
print(f"Final Energy: {final_energy}")

write('slab_opt.xyz', cu_slab)

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')

