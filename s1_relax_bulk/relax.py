# [1] https://doi.org/10.1063/1.4994149

from math import log
import pickle
from ase.build import bulk
from ase.io import write
from ase import Atoms
from ase.optimize import LBFGS
from ase.calculators.espresso import Espresso
from time import time


start = time()
logfile = 'ase.log'

initial_output_file = 'bulk_init.xyz'
final_output_file = 'bulk_final.xyz'
cu_bulk = bulk('Cu', crystalstructure='fcc', a=3.6, cubic=True)
write(initial_output_file, cu_bulk)

pseudopotentials = {'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF'}
input_settings = {
    'control': {
        'calculation': 'vc-relax',
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
    'cell': {
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'press_conv_thr': 0.5,
    }
}

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(4, 4, 4),
                pseudo_dir='/home/harris.se/espresso/pseudos/',
                input_data=input_settings)
cu_bulk.calc = calc
opt = LBFGS(cu_bulk, logfile=logfile, trajectory='qn.traj')
opt.run(fmax=0.001)

write(final_output_file, cu_bulk)

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')


