import os

import numpy as np
import matplotlib.pyplot as plt

import adlib.bulk.convergence
import adlib.bulk.eos


kpts_dir = '/work/westgroup/harris.se/espresso/adlib/test/bulk/kpts_converge'
kpts = [1, 2, 3, 4, 5, 6, 7]
lattice_constants = np.zeros(len(kpts))
for i, k in enumerate(kpts):
    calc_dir = os.path.join(kpts_dir, str(k))
    lattice_constants[i] = adlib.bulk.eos.plot_eos(calc_dir)


best_guess = lattice_constants[-1]
epsilon = 0.01
lower_lim = best_guess - epsilon
upper_lim = best_guess + epsilon

kpts_filename = os.path.join(kpts_dir, 'kpts_converge.png')
plt.plot(kpts, lattice_constants)
plt.xlabel('kpts')
plt.ylabel('lattice constant')

plt.axhline(y=lower_lim, color='r')
plt.axhline(y=upper_lim, color='r')
plt.savefig(kpts_filename)
