from ase.io.espresso import read_espresso_out
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.build import add_adsorbate
from matplotlib import pyplot as plt
import numpy as np
from time import time


start = time()
logfile = 'ase.log'
height = 1.0

# read in the results from the slab + vacuum relaxation
slab_file = '/home/harris.se/dft_adsorption/s2_relax2_slab_7p5/espresso.pwo'

with open(slab_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    clean_slab = traj[-1]

initial_energy = clean_slab.get_potential_energy()
print(f"Initial Energy: {initial_energy}")


# import optimized co2
co2_file = '/home/harris.se/dft_adsorption/s3_co2_7p5/espresso.pwo'
with open(co2_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    co2 = traj[-1]


# Calculator Settings
pseudopotentials = {'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
                    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
                    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'}
input_settings = {
    'control': {
        'calculation': 'scf',
        'forc_conv_thr': 0.001
    },
    'system': {
        'occupations': 'smearing',
        'degauss': 0.1,
        'ecutwfc': 50,
        'ecutrho': 500
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    },
    'cell': {
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'press_conv_thr': 0.5,
    }
}

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(4, 4, 1),
                pseudo_dir='/home/harris.se/espresso/pseudos/',
                input_data=input_settings)
clean_slab.calc = calc

# place the co2
co2_index = [atom.index for atom in co2 if atom.symbol == 'C']
add_adsorbate(clean_slab, co2, height=height, mol_index=co2_index[0])

write(f'system.xyz', clean_slab)
system_energy = clean_slab.get_potential_energy()
print(f'Height:\t{height}')
print(f'System Energy:\t{system_energy}')

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')

