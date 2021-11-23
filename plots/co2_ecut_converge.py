import numpy as np
from matplotlib import pyplot as plt


energies = [
    -1369.3774798674665,
    -1384.2167264170578,
    -1384.635030422279,
    -1384.6503690711932
]

ecuts = [
    25,
    36,
    47,
    58
]

runtimes = np.array([
    10811.015125513077,
    29409.5446434021,
    23436.005511045456,
    32461.756583929062
])


plt.plot(ecuts, energies, marker='o')
plt.ylabel('Energy (eV)')
plt.xlabel('ecutwfc (Ry)')
plt.title('Convergence of Plane Wave Cutoff Energy for CO2 Geometry Optimization')
plt.show()

plt.scatter(ecuts, runtimes / 3600.0, color='C1')
plt.title('CO2 Geometry Optimization Runtime vs. Plane Wave Cutoff')
plt.xlabel('ecutwfc (Ry)')
plt.ylabel('Run time (hours)')
plt.show()
