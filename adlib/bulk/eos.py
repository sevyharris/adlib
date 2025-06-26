"""
Module for setting up and running bulk energy calculations in Quantum Espresso

produces equation of state - Energy vs. lattice constant

Expects/creates the following directory structure:

bulk
____bulk.pwo
____vc_relax
________relax_bulk.py
________espresso.pwo
________run.sh
____eos_coarse
________run_qe_jobs.sh
________run0
____________calc.py
____________espresso.pwo
________run1
____________calc.py
____________espresso.pwo
________runN
____________calc.py
____________espresso.pwo
____eos_fine
________run_qe_jobs.sh
________run0
____________calc.py
____________espresso.pwo
________run1
____________calc.py
____________espresso.pwo
________runN
____________calc.py
____________espresso.pwo
____convergence_check
________kpts
________ecutwfc
________smear


"""

import os
import sys
import glob

import numpy as np
from matplotlib import pyplot as plt
from ase.io.espresso import read_espresso_out
import ase.eos

import adlib.bulk.calc


def setup_eos(bulk_dir, lattice_constant_guess, metal='Cu', N=21, half_range=0.05, magnetism=None):
    """
    script to set up N jobs

    N = number of steps/jobs at which to compute enregy
    half_range is half of the domain of the lattice constant guess
        for example, if my guess for lattice constant is 4.0, and I want to
        compute energies for lattice constants ranging from 3.5 - 4.5,
        then half_range should be set to 0.5 Angstroms
    """
    eos_dir = os.path.join(bulk_dir, 'eos')
    os.makedirs(eos_dir, exist_ok=True)

    deltas = np.linspace(-half_range, half_range, N)
    lattice_constants = deltas + lattice_constant_guess

    for i, lattice_constant in enumerate(lattice_constants):
        calc_dir = os.path.join(eos_dir, f'run_{i:04}')
        adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, magnetism=magnetism, ecutwfc=1000, kpt=9, smear=0.1, nproc=16)

    adlib.bulk.calc.make_scf_run_file_array(eos_dir, i, job_name='bulk_eos')


def get_lattice_constant_from_atoms(atoms):
    fcc_metals = ['Pt', 'Ni']
    bcc_metals = ['Cr', 'Fe']

    r = atoms.get_distances(0, 1)[0] / 2.0
    if atoms[0].symbol in fcc_metals:
        lattice_constant = 2.0 * r * np.sqrt(2)
    elif atoms[0].symbol in bcc_metals:
        lattice_constant = 4.0 * r / np.sqrt(3.0)
    else:
        raise NotImplemented(f'Metal {atoms[0].symbol} not in fcc or bcc metals')
    return lattice_constant


def setup_eos_coarse(bulk_dir, lattice_constant_guess, metal='Cu', magnetism=None):
    """
    script to set up coarse calculation of E vs. lattice constant
    """
    eos_dir = os.path.join(bulk_dir, 'eos_coarse')
    os.makedirs(eos_dir, exist_ok=True)

    deltas = np.linspace(-0.05, 0.05, 11)
    lattice_constants = deltas + lattice_constant_guess

    for i, lattice_constant in enumerate(lattice_constants):
        calc_dir = os.path.join(eos_dir, f'run_{i:04}')
        adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, magnetism=magnetism)

    adlib.bulk.calc.make_scf_run_file_array(eos_dir, i, job_name='bulk_eos_coarse')


def setup_eos_fine(bulk_dir, lattice_constant_guess, metal='Cu', crystal_structure='fcc', magnetism=None):
    """
    script to set up fine calculation of E vs. lattice constant
    """
    eos_dir = os.path.join(bulk_dir, 'eos_fine')
    os.makedirs(eos_dir, exist_ok=True)

    deltas = np.linspace(-0.005, 0.005, 11)
    lattice_constants = deltas + lattice_constant_guess

    for i, lattice_constant in enumerate(lattice_constants):
        calc_dir = os.path.join(eos_dir, f'run_{i:04}')
        adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, crystal_structure=crystal_structure, magnetism=magnetism)

    adlib.bulk.calc.make_scf_run_file_array(eos_dir, i, job_name='bulk_eos_fine')


def run_eos(calc_dir):
    # expecting the eos directory
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(calc_dir)
    eos_job = job_manager.SlurmJob()
    cmd = "sbatch run_qe_jobs.sh"
    eos_job.submit(cmd)
    os.chdir(cur_dir)


def analyze_eos(calc_dir):
    """function to find the lowest energy lattice constant
    in a folder full of QE runs

    returns the minimum energy lattice constant using ase's eos fit
    assumes fcc geometry
    """

    pwo_files = glob.glob(os.path.join(calc_dir, '*', 'espresso.pwo'))
    N = len(pwo_files)
    pwo_files.sort()

    energies = []
    volumes = []

    for pwo_file in pwo_files:
        with open(pwo_file, 'r') as f:
            traj = list(read_espresso_out(f, index=slice(None)))
            try:
                atoms = traj[-1]
            except IndexError:
                continue
            energies.append(atoms.get_potential_energy())
            volumes.append(atoms.cell.volume)

    bulk_eos = ase.eos.EquationOfState(volumes, energies)
    v0, e0, B = bulk_eos.fit()
    a0 = np.float_power(v0, 1.0 / 3.0)  # change this to match bcc or fcc
    return a0


def plot_energy_vs_lattice(calc_dir, dest_dir=None):
    """function to plot energy vs. lattice constant
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
            try:
                atoms = traj[-1]
                energies.append(atoms.get_potential_energy())
            except IndexError:
                energies.append(np.nan)

            lattice_constant = get_lattice_constant_from_atoms(atoms)
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
    plt.savefig(os.path.join(dest_dir, 'energy_vs_lattice.png'))


def plot_eos(calc_dir, dest_dir=None):
    """function to plot energy vs. lattice constant, using ase EOS fitting
    """
    if dest_dir is None:
        dest_dir = calc_dir

    pwo_files = glob.glob(os.path.join(calc_dir, '*', 'espresso.pwo'))
    N = len(pwo_files)
    pwo_files.sort()

    energies = []
    lattice_constants = []
    volumes = []

    for pwo_file in pwo_files:
        with open(pwo_file, 'r') as f:
            try:
                traj = list(read_espresso_out(f, index=slice(None)))
                atoms = traj[-1]
            except IndexError:
                continue
            energies.append(atoms.get_potential_energy())
            lattice_constant = get_lattice_constant_from_atoms(atoms)
            lattice_constants.append(lattice_constant)
            volumes.append(atoms.cell.volume)

    bulk_eos = ase.eos.EquationOfState(volumes, energies)
    v0, e0, B = bulk_eos.fit()
    a0 = np.float_power(v0, 1.0 / 3.0)
    plot_fname = os.path.join(dest_dir, 'ase_eos.png')
    ax = bulk_eos.plot(filename=plot_fname, show=False)
    ax.clear()
    return a0
