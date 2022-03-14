"""
Module for checking bulk convergence vs. kpts, ecutwfc, and smearing
"""

import os
import sys
import glob

import numpy as np
from matplotlib import pyplot as plt
from ase.io.espresso import read_espresso_out
import adlib.bulk.calc


def setup_kpts_converge(bulk_dir, lattice_constant, metal='Cu'):
    """
    script to set up N jobs to check kpts convergence
    """
    # TODO check that ecutwfc and smear are reasonable values
    kpts_dir = os.path.join(bulk_dir, 'kpts_converge')
    os.makedirs(kpts_dir, exist_ok=True)

    kpts = [k for k in range(1, 22)]  # 1-21

    for i, k in enumerate(kpts):
        calc_dir = os.path.join(kpts_dir, f'run_{i:04}')
        adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, ecutwfc=1000, kpt=k, smear=0.1, nproc=16)

    adlib.bulk.calc.make_scf_run_file_array(kpts_dir, i, job_name='kpts_bulk_converge')


def setup_ecutwfc_converge(bulk_dir, lattice_constant, metal='Cu'):
    """
    script to set up N jobs to check ecutwfc convergence
    """
    # TODO check that kpts and smear are reasonable values
    ecutwfc_dir = os.path.join(bulk_dir, 'ecutwfc_converge')
    os.makedirs(ecutwfc_dir, exist_ok=True)

    ecuts = [40.0, 60.0, 80.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0, 1200.0, 1400.0, 1600.0, 1800.0, 2000.0, 3000.0]

    for i, ecut in enumerate(ecuts):
        calc_dir = os.path.join(ecutwfc_dir, f'run_{i:04}')
        adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, ecutwfc=ecut, kpt=7, smear=0.1, nproc=16)

    adlib.bulk.calc.make_scf_run_file_array(ecutwfc_dir, i, job_name='ecut_bulk_converge')


def setup_smear_converge(bulk_dir, lattice_constant, metal='Cu'):
    """
    script to set up N jobs to check MV smearing convergence
    """
    # TODO check that ecutwfc and kpts are reasonable values
    smear_dir = os.path.join(bulk_dir, 'smear_converge')
    os.makedirs(smear_dir, exist_ok=True)

    smears = [0.01, 0.025, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5]
    for i, smear in enumerate(smears):
        calc_dir = os.path.join(smear_dir, f'run_{i:04}')
        adlib.bulk.calc.make_scf_calc_file(calc_dir, lattice_constant, metal=metal, ecutwfc=1000, kpt=7, smear=smear, nproc=16)

    adlib.bulk.calc.make_scf_run_file_array(smear_dir, i, job_name='smear_bulk_converge')
