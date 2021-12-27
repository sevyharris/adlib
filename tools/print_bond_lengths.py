import sys
from ase.io.espresso import read_espresso_out
from ase import Atoms

pwo_file = sys.argv[1]

with open(pwo_file, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    atoms = traj[-1]
    for i in range(1, len(atoms)):
        print(atoms.get_distances(0, i))
    print(atoms.get_potential_energy())

