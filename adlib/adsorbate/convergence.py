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
from ase.io import read, write
import adlib.adsorbate.calc


def setup_relax_adsorbate(adsorbate_dir, xyz_dir=None):
    ads_name = os.path.basename(adsorbate_dir)
    if xyz_dir is None:
        # TODO make this relative to the package and ship the code with some example adsorbate xyzs
        xyz_dir = '/work/westgroup/harris.se/espresso/qe_workflow/resources/adsorbates/'
    os.makedirs(adsorbate_dir, exist_ok=True)
    xyz_file = os.path.join(xyz_dir, f'{ads_name}.xyz')
    shutil.copy(xyz_file, adsorbate_dir)
    make_relax_ads_script(adsorbate_dir)
    make_run_relax_ads_script(adsorbate_dir)


def setup_vacuum_converge(adsorbate_dir, xyz_dir=None, adsorbate_pwo=None):
    """
    script to set up N jobs to check vacuum convergence
    """
    ads_name = os.path.basename(adsorbate_dir)
    vacuum_dir = os.path.join(adsorbate_dir, 'vacuum_converge')
    os.makedirs(vacuum_dir, exist_ok=True)

    if adsorbate_pwo:
        with open(adsorbate_pwo, 'r') as f:
            traj = list(read_espresso_out(f, index=slice(None)))
            adsorbate = traj[-1]
    elif xyz_dir:
        xyz_file = os.path.join(xyz_dir, f'{ads_name}.xyz')
        adsorbate = read(xyz_file)
    else:
        # TODO make this relative to the package and ship the code with some example adsorbate xyzs
        xyz_dir = '/work/westgroup/harris.se/espresso/qe_workflow/resources/adsorbates/'
        xyz_file = os.path.join(xyz_dir, f'{ads_name}.xyz')
        adsorbate = read(xyz_file)

    vacuums = [v for v in range(5, 16)]

    for i, vac in enumerate(vacuums):
        calc_dir = os.path.join(vacuum_dir, f'run_{i:04}')
        ads_filename = os.path.join(calc_dir, ads_name + '.xyz')
        write(ads_filename, adsorbate)
        adlib.adsorbate.calc.make_relax_ads_script(calc_dir, vacuum=vac)

    adlib.adsorbate.calc.make_run_array(vacuum_dir, i, job_name='ads_vacuum')


def run_vacuum_converge(adsorbate_dir):
    import job_manager
    cur_dir = os.getcwd()
    vacuum_dir = os.path.join(adsorbate_dir, 'vacuum_converge')
    os.chdir(vacuum_dir)
    ads_job = job_manager.SlurmJob()
    cmd = "sbatch run_qe_jobs.sh"
    ads_job.submit(cmd)
    os.chdir(cur_dir)
