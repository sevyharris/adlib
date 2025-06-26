"""
Module for checking bulk convergence vs. kpts, ecutwfc, and smearing
____convergence_check
________kpts
____________1
________________run_0000
____________________calc.py
____________________espresso.pwo
________________run_0001
____________________calc.py
____________________espresso.pwo
________________run_000N
____________________calc.py
____________________espresso.pwo
____________2
________________run_0000
____________________calc.py
____________________espresso.pwo
________________run_0001
____________________calc.py
____________________espresso.pwo
________________run_000N
____________________calc.py
____________________espresso.pwo
____________N
________________run_0000
____________________calc.py
____________________espresso.pwo
________________run_0001
____________________calc.py
____________________espresso.pwo
________________run_000N
____________________calc.py
____________________espresso.pwo
________ecutwfc
________smear
"""

import os
import sys
import glob

import numpy as np
from matplotlib import pyplot as plt
from ase.io.espresso import read_espresso_out
import adlib.bulk.calc


def setup_kpts_converge(bulk_dir, lattice_constant_center, metal='Cu', crystal_structure='fcc'):
    """
    script to set up N jobs to check kpts convergence
    """

    # TODO check that ecutwfc and smear are reasonable values
    kpts_dir = os.path.join(bulk_dir, 'kpts_converge')
    os.makedirs(kpts_dir, exist_ok=True)

    kpts = [k for k in range(1, 22)]  # 1-21
    lattice_constants = lattice_constant_center + np.linspace(-0.05, 0.05, 11)

    for i, k in enumerate(kpts):
        for j, lattice_constant in enumerate(lattice_constants):
            calc_dir = os.path.join(kpts_dir, str(k), f'run_{j:04}')
            adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, crystal_structure=crystal_structure, ecutwfc=500, kpt=k, smear=0.1, nproc=16)

        adlib.bulk.calc.make_scf_run_file_array(os.path.join(kpts_dir, str(k)), j, job_name='kpts_bulk_converge')


def setup_ecutwfc_converge(bulk_dir, lattice_constant_center, metal='Cu', crystal_structure='fcc'):
    """
    script to set up N jobs to check ecutwfc convergence
    """
    # TODO check that kpts and smear are reasonable values
    ecutwfc_dir = os.path.join(bulk_dir, 'ecutwfc_converge')
    os.makedirs(ecutwfc_dir, exist_ok=True)

    lattice_constants = lattice_constant_center + np.linspace(-0.05, 0.05, 11)
    ecuts = [40.0, 60.0, 80.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 3000.0]

    for i, ecut in enumerate(ecuts):
        for j, lattice_constant in enumerate(lattice_constants):
            calc_dir = os.path.join(ecutwfc_dir, str(ecut), f'run_{j:04}')
            adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, crystal_structure=crystal_structure, ecutwfc=ecut, kpt=7, smear=0.1, nproc=16)

        adlib.bulk.calc.make_scf_run_file_array(os.path.join(ecutwfc_dir, str(ecut)), j, job_name='ecut_bulk_converge')


def setup_smear_converge(bulk_dir, lattice_constant_center, metal='Cu', crystal_structure='fcc'):
    """
    script to set up N jobs to check MV smearing convergence
    """
    # TODO check that ecutwfc and kpts are reasonable values
    smear_dir = os.path.join(bulk_dir, 'smear_converge')
    os.makedirs(smear_dir, exist_ok=True)

    smears = [0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5]
    lattice_constants = lattice_constant_center + np.linspace(-0.05, 0.05, 11)
    for i, smear in enumerate(smears):
        for j, lattice_constant in enumerate(lattice_constants):
            calc_dir = os.path.join(smear_dir, str(smear), f'run_{j:04}')
            adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, crystal_structure=crystal_structure, ecutwfc=1000, kpt=7, smear=smear, nproc=16)

        adlib.bulk.calc.make_scf_run_file_array(os.path.join(smear_dir, str(smear)), j, job_name='smear_bulk_converge')


def plot_kpts_converge(calc_dir, dest_dir=None):
    """function to plot lattice constant vs. kpts
    """
    if dest_dir is None:
        dest_dir = calc_dir

    pwo_files = glob.glob(os.path.join(calc_dir, '*', 'espresso.pwo'))
    N = len(pwo_files)
    pwo_files.sort()

    energies = []
    lattice_constants = []

    for pwo_file in pwo_files:
        with open(pwo_file, 'r') as f:
            traj = list(read_espresso_out(f, index=slice(None)))
            atoms = traj[-1]
            energies.append(atoms.get_potential_energy())
            lattice_constant = atoms.get_distances(0, 1)[0] * np.sqrt(2)
            lattice_constants.append(lattice_constant)

    fig, ax = plt.subplots()
    plt.plot(lattice_constants, energies, marker='o')

    # label the minimum
    label_min = True
    if label_min:
        min_energy = np.min(energies)
        min_i = energies.index(min_energy)
        ax.annotate(
            f'({np.round(lattice_constants[min_i], 3)}, {np.round(min_energy, 3)})',
            xy=(lattice_constants[min_i], min_energy),
            xytext=(lattice_constants[min_i], np.mean(energies)),
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
        )

    plt.ylabel('Energy (eV)')
    ax.yaxis.get_major_formatter().set_useOffset(False)
    plt.xlabel(r'Lattice Constant ($\AA$)')
    plt.title('Bulk Energy vs. Lattice Constant')
    plt.savefig(os.path.join(dest_dir, 'equation_of_state.png'))


def plot_ecutwfc_converge(calc_dir):
    pass


def plot_smear_converge(calc_dir):
    pass
