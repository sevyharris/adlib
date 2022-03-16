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
    lattice_constants[i] = adlib.bulk.eos.plot_eos2(calc_dir)

kpts_filename = os.path.join(kpts_dir, 'kpts_converge.png')
plt.plot(kpts, lattice_constants)
plt.savefig(kpts_filename)
