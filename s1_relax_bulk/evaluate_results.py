from ase.io.espresso import read_espresso_out
from ase import Atoms
from ase.calculators.espresso import Espresso
from ase.io.espresso import write_espresso_in


with open('espresso.pwo', 'r') as f:
    traj = list(read_espresso_out(f, index=slice(None)))
    print(traj)
#my_structure = list(read_espresso_out('espresso.pwo'))
#for atom in my_structure:
#    print(atom)

for structure in traj:
    print(structure.get_volume())
#print(my_structure())
#print(my_structure.get_volume())
#print(my_structure.cell.get_bravais_lattice())

