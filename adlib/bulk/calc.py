"""
Module for setting up and running bulk energy calculations in Quantum Espresso

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



"""

import os
import adlib.env


environment = adlib.env.load_environment()


def make_scf_run_file(calc_dir, nproc=16, job_name='bulk_energy'):
    bash_filename = os.path.join(calc_dir, 'run.sh')
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')

        if environment == 'DISCOVERY':
            f.write('#SBATCH --time=24:00:00\n')
            f.write(f'#SBATCH --job-name={job_name}\n')
            f.write('#SBATCH --mem=40Gb\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --ntasks={nproc}\n')
            f.write('#SBATCH --partition=short,west\n')
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')
        elif environment == 'EXPLORER':
            f.write('#SBATCH --time=24:00:00\n')
            f.write(f'#SBATCH --job-name={job_name}\n')
            f.write('#SBATCH --mem=40Gb\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --ntasks={nproc}\n')
            f.write('#SBATCH --partition=short,west\n\n')
            f.write('module load OpenMPI/4.1.6\n')

        # f.write(f'cd {calc_dir}\n')
        f.write(f'python calc.py\n')


def make_scf_run_file_array(dest_dir, N_runs, job_name='bulk_energy', nproc=16):
    bash_filename = os.path.join(dest_dir, 'run_qe_jobs.sh')
    run_i_dir = os.path.abspath(os.path.join(dest_dir, 'run_$RUN_i'))
    # write the array job file
    with open(bash_filename, 'w') as f:
        # TODO add array option for single node environment
        f.write('#!/bin/bash\n\n')
        f.write('#SBATCH --time=24:00:00\n')
        f.write(f'#SBATCH --job-name={job_name}\n')
        f.write('#SBATCH --mem=40Gb\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write(f'#SBATCH --ntasks={nproc}\n')
        f.write('#SBATCH --partition=short,west\n')
        f.write(f'#SBATCH --array=0-{N_runs}\n\n')

        f.write('# create a variable that includes leading zeros\n')
        f.write('RUN_i=$(printf "%04.0f" $SLURM_ARRAY_TASK_ID)\n')

        if environment == 'DISCOVERY':
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')
        elif environment == 'EXPLORER':
            f.write('module load OpenMPI/4.1.6\n')

        f.write(f'cd {run_i_dir}\n')
        f.write(f'python calc.py\n')


def make_scf_calc_file(calc_dir, lattice_constant, metal='Cu', crystal_structure='fcc', ecutwfc=1000, ecutrho=None, kpt=9, smear=0.1, nproc=16, magnetism=None, mixing_beta=0.7):

    magnetism_line = ""
    if magnetism is not None:
        magnetism_line = f'initial_magnetic_moments = {magnetism}' + \
            '\nassert len(bulk) == len(initial_magnetic_moments)' + \
            '\nbulk.set_initial_magnetic_moments(initial_magnetic_moments)'
    ecutrho_line = ""
    if ecutrho is not None:
        ecutrho_line = f"        'ecutrho': {ecutrho},"

    electrons_line = ""
    if mixing_beta != 0.7:
        electrons_line = """
    'electrons': {
        'mixing_beta': """ + str(mixing_beta) + """,
        'electron_maxstep': 200,
    },"""

    python_file_lines = [
        "import os",
        "import sys",
        "import time",
        "from ase.calculators.espresso import Espresso, EspressoProfile",
        "import ase.build",
        "",
        "",
        "T = time.localtime()",
        "start = time.time()",
        "logfile = 'ase.log'",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Start: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "",
        "",
        "espresso_settings = {",
        "    'control': {",
        "        'verbosity': 'high',",
        "        'calculation': 'scf',",
        "        'disk_io': 'none',",
        "    },",
        "    'system': {",
        "        'input_dft': 'BEEF-VDW',",
        "        'occupations': 'smearing',",
        "        'smearing': 'mv',",
        f"        'degauss': {smear},",
        f"        'ecutwfc': {ecutwfc},",
        ecutrho_line,
        "    },",
        electrons_line,
        "}",
        "",
        "",
        f"bulk = ase.build.bulk('{metal}', crystalstructure='{crystal_structure}', a={lattice_constant}, cubic=True)",
        "",
        magnetism_line,
        "pw_executable = os.environ['PW_EXECUTABLE']",
        "",
        "pseudopotentials = {",
        "    'C': 'C_ONCV_PBE-1.2.upf',",
        "    'Cu': 'Cu_ONCV_PBE-1.2.upf',",
        "    'O': 'O_ONCV_PBE-1.2.upf',",
        "    'N': 'N_ONCV_PBE-1.2.upf',",
        "    'H': 'H_ONCV_PBE-1.2.upf',",
        "    'Pt': 'Pt_ONCV_PBE-1.2.upf',",
        "    'Pd': 'Pd_ONCV_PBE-1.2.upf',",
        "    'Au': 'Au_ONCV_PBE-1.2.upf',",
        "    'Ag': 'Ag_ONCV_PBE-1.2.upf',",
        "    'Al': 'Al_ONCV_PBE-1.2.upf',",
        "    'Ni': 'Ni_ONCV_PBE-1.2.upf',",
        "    'Fe': 'Fe_ONCV_PBE-1.2.upf',",
        "    'Cr': 'Cr_ONCV_PBE-1.2.upf',",
        "}",
        "",
        f"command = f'mpirun -np {nproc} " + "{pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "print(command)",
        "profile = EspressoProfile(argv=command.split())",
        "",
        "espresso = Espresso(",
        "    # command=command,",
        "    profile=profile,",
        "    pseudopotentials=pseudopotentials,",
        "    tstress=True,",
        "    tprnfor=True,",
        f"    kpts=({kpt}, {kpt}, {kpt}),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "bulk.calc = espresso",
        "energy = bulk.get_potential_energy()",
        "",
        "",
        "T = time.localtime()",
        "end = time.time()",
        "duration = end - start",
        "",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Energy: {energy} eV\\n')",
        "    f.write(f'End: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
        "",
    ]
    os.makedirs(calc_dir, exist_ok=True)
    calc_filename = os.path.join(calc_dir, f'calc.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])
