"""
Module for setting up and running vc-relax in Quantum Espresso

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
import adlib.env


environment = adlib.env.load_environment()


def setup_vc_relax(bulk_dir, metal='Cu', lattice_constant_guess=3.6, crystal_structure='fcc', nproc=16):
    os.makedirs(bulk_dir, exist_ok=True)
    vc_relax_dir = os.path.join(bulk_dir, 'vc_relax')
    make_vc_relax_script(vc_relax_dir, metal, lattice_constant_guess, crystal_structure=crystal_structure, nproc=nproc)
    make_run_vc_relax_script(vc_relax_dir)


def make_run_vc_relax_script(calc_dir, nproc=16):
    """
    Make a shell script to run the python script
    """
    os.makedirs(calc_dir, exist_ok=True)
    bash_filename = os.path.join(calc_dir, 'run.sh')
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        if environment == 'DISCOVERY':
            f.write('#SBATCH --time=24:00:00\n')
            f.write('#SBATCH --job-name=vc_relax\n')
            f.write('#SBATCH --mem=40Gb\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --ntasks={nproc}\n')
            f.write('#SBATCH --partition=short,west\n')
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')
        # f.write(f'cd {calc_dir}\n')
        f.write(f'python relax_bulk.py\n')


def make_vc_relax_script(calc_dir, metal, lattice_constant, crystal_structure='fcc', magnetism=None, nproc=16):
    """
    Make a python script that uses ase to run quantum espresso
    """

    python_file_name = os.path.join(calc_dir, 'relax_bulk.py')
    kpts = 5
    ecutwfc = 500

    magnetism_line = ""
    if str(magnetism).lower() == 'antiferromagnetic':
        magnetism_line = "bulk.set_initial_magnetic_moments([1.0, -1.0])\nassert len(bulk) == 2"
    elif str(magnetism).lower() == 'ferromagnetic':
        magnetism_line = "bulk.set_initial_magnetic_moments([1.0, 1.0])\nassert len(bulk) == 2"

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
        "        'calculation': 'vc-relax',",
        "    },",
        "    'system': {",
        "        'input_dft': 'BEEF-VDW',",
        "        'occupations': 'smearing',",
        "        'smearing': 'mv',",
        "        'degauss': 0.1,",
        f"        'ecutwfc': {ecutwfc},",
        "    },",
        "    'ions': {",
        "        'ion_dynamics': 'bfgs',",
        "    },",
        "    'cell': {",
        "        'cell_dynamics': 'bfgs',",
        "        'press': 0.0,",
        "        'press_conv_thr': 0.005,",
        "    }",
        "}",
        "",
        "",
        f"bulk_metal = ase.build.bulk('{metal}', crystalstructure='{crystal_structure}', a={lattice_constant}, cubic=True)",
        magnetism_line,
        "",
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
        f"    kpts=({kpts}, {kpts}, {kpts}),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "bulk_metal.calc = espresso",
        "energy = bulk_metal.get_potential_energy()",
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
    with open(python_file_name, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])


def run_vc_relax(bulk_dir):
    import job_manager
    calc_dir = os.path.join(bulk_dir, 'vc_relax')
    cur_dir = os.getcwd()
    os.chdir(calc_dir)
    if environment in ['DISCOVERY', 'EXPLORER']:
        vc_relax_job = job_manager.SlurmJob()
        cmd = "sbatch run.sh"
    elif environment == 'SINGLE_NODE':
        vc_relax_job = job_manager.DefaultJob()
        cmd = "/bin/bash run.sh"
    elif environment == 'THETA':
        raise NotImplementedError
    else:
        raise NotImplementedError

    vc_relax_job.submit(cmd)
    os.chdir(cur_dir)
