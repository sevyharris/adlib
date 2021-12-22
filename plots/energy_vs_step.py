import numpy as np
from matplotlib import pyplot as plt
from ase.io.espresso import read_espresso_out
from ase import Atoms

pwo_file1 = '/home/moon/dft_adsorption/plots/slab1.pwo'
pwo_file2 = '/home/moon/dft_adsorption/plots/slab2.pwo'

energies = []
# lattice_constants = []

with open(pwo_file1, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))

for atoms in traj:
    # lattice_constants.append(atoms.cell[0][0])
    energies.append(atoms.get_potential_energy())

with open(pwo_file2, 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))

for atoms in traj:
    # lattice_constants.append(atoms.cell[0][0])
    energies.append(atoms.get_potential_energy())

fig, ax = plt.subplots()
plt.plot(energies, marker='o')
plt.ylabel('Energy (eV)')
ax.yaxis.get_major_formatter().set_useOffset(False)
plt.xlabel(r'Lattice Constant ($\AA$)')
plt.title('Bulk Energy vs. Lattice Constant')
plt.show()

print(energies)
# print(lattice_constants)
