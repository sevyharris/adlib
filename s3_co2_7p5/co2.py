from ase.calculators.espresso import Espresso
from ase.optimize import LBFGS
from ase import Atoms
from ase.io import write
from time import time


start = time()
logfile = 'ase.log'


co2_init_file = 'co2_init.xyz'
co2_opt_file = 'co2_opt.xyz'

pseudopotentials = {'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
                    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'}

input_settings = {
    'control': {
        'calculation': 'relax',
        'forc_conv_thr': 0.001
    },
    'system': {
        'ecutwfc': 50,  # energy cutoffs from [1]
        'ecutrho': 500
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    },
}

d = 1.0
co2 = Atoms('CO2', positions=[(0, 0, 0), (d, 0, 0), ((-d, 0, 0))])
co2.center(vacuum=7.5)
write(co2_init_file, co2)

calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=None,
                pseudo_dir='/home/harris.se/espresso/pseudos/',
                input_data=input_settings)
co2.calc = calc

opt = LBFGS(co2, trajectory='co2.traj', logfile=logfile)
opt.run(fmax=0.005)

write(co2_opt_file, co2)

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')

