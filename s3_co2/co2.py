from ase.calculators.espresso import Espresso
from ase.optimize import LBFGS
from ase import Atoms
from time import time


start = time()
logfile = 'ase.log'

pseudopotentials = {'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
                    'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF'}

d = 2.0
co2 = Atoms('CO2', positions=[(0, 0, 0), (0, 0, d), ((0, 0, -d))])
co2.center(vacuum=12)
calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(3, 3, 3),
                pseudo_dir='/home/harris.se/espresso/pseudos/')
co2.calc = calc

opt = LBFGS(co2, trajectory='co2.traj', logfile=logfile)
opt.run(fmax=0.005)

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds')

