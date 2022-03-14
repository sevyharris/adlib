import sys
from ase.io.espresso import read_espresso_out
from ase.io.trajectory import Trajectory
from ase import Atoms
import os

input_file = sys.argv[1]

if not os.path.exists(input_file):
    raise OSError

if input_file.endswith('.pwo'):
    with open(input_file, 'r') as f:
        traj = list(read_espresso_out(f, index=slice(None)))
elif input_file.endswith('.traj'):
    traj = Trajectory(input_file)

else:
    raise ValueError('Input file type not yet implemented')


energy = traj[-1].get_potential_energy()
print(energy)
                                                                                  
