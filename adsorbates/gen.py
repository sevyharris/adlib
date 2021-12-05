from ase.io import write
from ase import Atoms


# CO2 - need to make sure the orientation is relative to the z normal of the slab
d = 1.0
co2 = Atoms('CO2', positions=[(0, 0, 0), (d, 0, 0), ((-d, 0, 0))])
write('CO2.xyz', co2)
