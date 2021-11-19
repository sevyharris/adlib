from ase.io.espresso import read_espresso_out
from ase import Atoms
from ase.io import read, write
from ase.calculators.espresso import Espresso
from ase.build import add_adsorbate
from matplotlib import pyplot as plt
import numpy as np
import copy


# import optimized slab
slab_file = 'slab.pwo'
with open(slab_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
# print(len(traj))
clean_slab = traj[-1]
write('opt_slab.xyz', clean_slab)

# import optimized co2
d = 2.0
co2 = Atoms('CO2', positions=[(0, 0, 0), (d, 0, 0), ((-d, 0, 0))])
co2.center(vacuum=12)
# co2_file = 'co2.pwo'
# with open(co2_file, 'r') as f:
#     traj = list(read_espresso_out(f, index=slice(None)))
# co2 = traj[-1]
# write('opt_co2.xyz', co2)

# place the co2
index = [atom.index for atom in co2 if atom.symbol == 'C']

heights = np.linspace(0.5, 5, 20)


pseudopotentials = {'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
                    'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
                    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'}
input_settings = {
    'control': {
        'calculation': 'scf',
        'forc_conv_thr': 0.001
    },
    'system': {
        'occupations': 'smearing',  # required for metals
        'degauss': 0.1,
        'ecutwfc': 36,  # energy cutoffs from [1]
        'ecutrho': 400
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
                tstress=True, tprnfor=True, kpts=(1, 1, 1),
                pseudo_dir='../pseudos/',
                input_data=input_settings)
clean_slab.calc = calc

for i, height in enumerate(heights):
    slab = copy.deepcopy(clean_slab)
    add_adsorbate(slab, co2, height=height, mol_index=index[0])
    write(f'system{i}.xyz', slab)
    print(f'h={height}', f'U={slab.get_potential_energy()}')
