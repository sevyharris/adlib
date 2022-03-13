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

# Script to set up energy and run an energy calculation in ase/quantum espresso
import os
import sys

import numpy as np



def setup_eos_coarse(bulk_dir, lattice_constant):
    pass

def setup_eos_fine(bulk_dir, lattice_constant):
    pass


def setup_energy_calc(calc_dir, lattice_constant, array_job=True):
    os.makedirs(calc_dir, exist_ok=True)
    make_scf_calc_file(calc_dir, lattice_constant, ecutwfc=1000, kpt=9, smear=0.1, nproc=16)
    if not array_job:
        make_scf_run_file(calc_dir)


def make_scf_run_file(calc_dir, nproc=16):
    bash_filename = os.path.join(calc_dir, 'run.sh')
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('#SBATCH --time=24:00:00\n')
        f.write('#SBATCH --job-name=bulk_energy\n')
        f.write('#SBATCH --mem=40Gb\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write(f'#SBATCH --ntasks={nproc}\n')
        f.write('#SBATCH --partition=short,west\n')
        f.write('module load gcc/10.1.0\n')
        f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
        f.write('module load scalapack/2.1.0-skylake\n\n')
        f.write(f'cd {calc_dir}\n')
        f.write(f'python calc.py\n')


def make_scf_run_file_array(dest_dir, N_runs):
    bash_filename = os.path.join(dest_dir, 'run_qe_jobs.sh')
    run_i_dir = os.path.abspath(os.path.join(dest_dir, 'run$SLURM_ARRAY_TASK_ID'))
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('#SBATCH --time=24:00:00\n')
        f.write('#SBATCH --job-name=bulk_energy\n')
        f.write('#SBATCH --mem=40Gb\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write('#SBATCH --ntasks=16\n')
        f.write('#SBATCH --partition=short,west\n')
        f.write(f'#SBATCH --array=0-{N_runs - 1}\n\n')
        f.write('module load gcc/10.1.0\n')
        f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
        f.write('module load scalapack/2.1.0-skylake\n\n')
        f.write(f'cd {run_i_dir}\n')
        f.write(f'python calc.py\n')

def make_scf_calc_file(calc_dir, lattice_constant, metal='Cu', ecutwfc=1000, kpt=9, smear=0.1, nproc=16):

    python_file_lines = [
        "import os",
        "import sys",
        "from time import time",
        "from ase.calculators.espresso import Espresso",
        "from ase.build import bulk",
        "",
        "",
        "start = time()",
        "logfile = 'ase.log'",
        "",
        "",
        "espresso_settings = {",
        "    'control': {",
        "        'verbosity': 'high',",
        "        'calculation': 'scf',",
        "    },",
        "    'system': {",
        "        'input_dft': 'BEEF-VDW',",
        "        'occupations': 'smearing',",
        "        'smearing': 'mv',",
        f"        'degauss': {smear},",
        f"        'ecutwfc': {ecutwfc},",
        "    },",
        "}",
        "",
        "",
        f"cu_bulk = bulk('{metal}', crystalstructure='fcc', a={lattice_constant}, cubic=True)",
        "",
        "pw_executable = os.environ['PW_EXECUTABLE']",
        "",
        "pseudopotentials = {",
        "    'C': 'C_ONCV_PBE-1.2.upf',",
        "    'Cu': 'Cu_ONCV_PBE-1.2.upf',",
        "    'O': 'O_ONCV_PBE-1.2.upf',",
        "    'N': 'N_ONCV_PBE-1.2.upf',",
        "    'H': 'H_ONCV_PBE-1.2.upf',",
        "}",
        "",
        f"command = f'mpirun -np {nproc} " + "{pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "print(command)",
        "",
        "espresso = Espresso(",
        "    command=command,",
        "    pseudopotentials=pseudopotentials,",
        "    tstress=True,",
        "    tprnfor=True,",
        f"    kpts=({kpt}, {kpt}, {kpt}),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "cu_bulk.calc = espresso",
        "energy = cu_bulk.get_potential_energy()",
        "",
        "",
        "end = time()",
        "duration = end - start",
        "",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Energy: {energy} eV\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
        "",
    ]

    calc_filename = os.path.join(calc_dir, f'calc.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])
