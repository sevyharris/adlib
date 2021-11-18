from ase.io.espresso import read_espresso_out
from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.io.espresso import write_espresso_in
from matplotlib import pyplot as plt
import numpy as np


with open('espresso.pwo', 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))

potential_energy = []
lattice_constants = []

for structure in traj:
    potential_energy.append(structure.get_potential_energy())
    lattice_constants.append(structure.cell[0][0])

print(potential_energy[-1])
exit(0)

a_start = np.round(lattice_constants[0], decimals=3)
steps = [x for x in range(0, len(traj))]
fig, ax = plt.subplots()
ax.plot(steps, potential_energy)
ax.yaxis.get_major_formatter().set_useOffset(False)
plt.title(f'Cu Bulk Energy, start a={a_start}')
plt.xlabel('Step')
plt.ylabel('Potential Energy')
plt.savefig('U.png')
plt.show()


plt.plot(steps, lattice_constants)
plt.title(f'Cu Bulk lattice constant, start a={a_start}')
plt.xlabel('Step')
plt.ylabel('Lattice constant (A)')
plt.savefig('a.png')
plt.show()
