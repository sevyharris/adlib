from ase.io.espresso import read_espresso_out
from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.io import write
from ase.io.trajectory import Trajectory

#traj = Trajectory('run2/co2.traj')
#atoms = traj[-1]
#write('co2.xyz', atoms)

#print(atoms)
#print(atoms.get_potential_energy())
#exit(0)

with open('espresso.pwo', 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))

potential_energy = []
lattice_constants = []

for structure in traj:
    potential_energy.append(structure.get_potential_energy())
    lattice_constants.append(structure.cell[0][0])

#print(lattice_constants[-1])
print(potential_energy[-1])
print(traj[-1].get_positions())
# print(traj[-1].get_distance(0, 1))

# write('co2.xyz', traj[-1])

