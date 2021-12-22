import numpy as np
from matplotlib import pyplot as plt


energies = [
    -1383.364832547786,     # 2
    -1384.1512250786136,    # 3
    -1384.2169148558912,    # 4
    -1384.210023845041,     # 5
    -1384.2151114214255,    # 6
    -1384.2189931253338,    # 7
    -1384.2190257789944,    # 8
    -1384.220851390738,     # 9
    -1384.2163017834127,    # 10
    -1384.2178952820518,    # 11
    -1384.2167771662887,    # 12
    -1384.218213519186,     # 13
    -1384.2183460386254,    # 14
    -1384.2174233005987,    # 15
]

vacuum = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]


runtimes = np.array([
    67.45378303527832,      # 2
    511.69579815864563,     # 3
    1430.513635635376,      # 4
    4907.019116401672,      # 5
    4559.0850439071655,     # 6
    7769.58563542366,       # 7
    15592.38325214386,      # 8
    16249.806287050247,     # 9
    22242.361535549164,     # 10
    5725.522533893585,      # 11
    31941.34457564354,      # 12
    10369.260548353195,     # 13
    21204.72437119484,      # 14
    52469.173446416855,     # 15
])


fig, ax = plt.subplots()
ax.yaxis.get_major_formatter().set_useOffset(False)
plt.plot(vacuum, energies, marker='o')
plt.ylabel('Energy (eV)')
plt.xlabel(r'Vacuum ($\AA$)')
plt.title('CO2 Geometry Optimization Final Energy vs. Vacuum Size')
plt.show()

# plt.plot(np.abs(np.diff(energies)), marker='o')
# plt.ylabel(r'$\Delta$Energy (eV)')
# plt.xlabel(r'Vacuum ($\AA$)')
# plt.title(r'CO2 Geometry Optimization $\Delta$ Energy vs. Vacuum Size')
# plt.show()

# print(np.diff(energies))

# plt.scatter(vacuum, runtimes / 3600.0, color='C1')
# plt.title('CO2 Geometry Optimization Runtime vs. Vacuum Size')
# plt.xlabel(r'Vacuum ($\AA$)')
# plt.ylabel('Run time (hours)')
# plt.show()
