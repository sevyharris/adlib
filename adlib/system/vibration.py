"""
Module for setting up and running frequency analysis in Quantum Espresso

system
____relax_system.py
____espresso.pwo
____vib
________espresso.pwo
"""

import os
import adlib.env


environment = adlib.env.load_environment()


def make_run_vib_analysis_script(calc_dir, nproc=16, job_name='vib'):
    bash_filename = os.path.join(calc_dir, 'run.sh')
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        if environment == 'DISCOVERY':
            # f.write('#SBATCH --time=24:00:00\n')
            f.write('#SBATCH --time=14-00:00:00\n')
            f.write(f'#SBATCH --job-name={job_name}\n')
            f.write('#SBATCH --mem=40Gb\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --ntasks={nproc}\n')
            f.write('#SBATCH --partition=west\n\n')
            # f.write('#SBATCH --partition=short\n')
            # f.write('#SBATCH --constraint=cascadelake\n\n')
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')

        f.write(f'python vib_analysis.py\n')


def make_vib_analysis_script(calc_dir, ecutwfc=50, kpt=5, smear=0.1, nproc=16, low_mixing_beta=False):
    """Function to make a python script to relax the slab-adsorption system
    """
    fmax = 0.01
    vacuum = 7.5

    electrons_str = ""
    if low_mixing_beta:
        electrons_str = """
    'electrons': {
        'mixing_beta': 0.3,
        'electron_maxstep': 200,
    },"""

    python_file_lines = [
        "import os",
        "import sys",
        "import time",
        "import numpy as np",
        "from ase.calculators.espresso import Espresso, EspressoProfile",
        "from ase.constraints import FixAtoms",
        "from ase.io.espresso import read_espresso_out",
        "from ase.vibrations import Vibrations",
        "",
        "",
        "T = time.localtime()",
        "start = time.time()",
        "",
        "logfile = 'ase.log'",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Start: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "",
        "# Read in the system espresso file",
        "with open('system.pwo', 'r') as f:",
        "    traj = list(read_espresso_out(f, index=slice(None)))",
        "    system = traj[-1]",
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
        "    },",
        f"{electrons_str}",
        "}",
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
        "",
        "# get the indicies of the gas molecules by assuming metal is most populous",
        "vals, counts = np.unique(system.get_atomic_numbers(), return_counts=True)",
        "metal_number = vals[np.argmax(counts)]",
        "adsorbate_indices = []",
        "for i, atom in enumerate(system):",
        "    if atom.number != metal_number:",
        "        adsorbate_indices.append(atom.index)",
        "",
        "pw_executable = os.environ['PW_EXECUTABLE']",
        f"command = f'mpirun -np {nproc} " + "{pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "profile = EspressoProfile(argv=command.split())",
        "espresso = Espresso(",
        "    # command=command,",
        "    profile=profile,",
        "    pseudopotentials=pseudopotentials,",
        "    tstress=True,",
        "    tprnfor=True,",
        f"    kpts=({kpt}, {kpt}, 1),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "system.calc = espresso",
        "",
        "vib = Vibrations(system, indices=adsorbate_indices)",
        "vib.run()",
        "vib.summary()",
        "freq = vib.get_frequencies()",
        "print(freq)",
        "print(vib.get_energies())",
        "",
        "T = time.localtime()",
        "end = time.time()",
        "duration = end - start",
        "",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Frequencies: {freq}')",
        "    f.write(f'End: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
    ]

    os.makedirs(calc_dir, exist_ok=True)
    calc_filename = os.path.join(calc_dir, f'vib_analysis.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])


def run_vib_analysis(vib_dir):
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(vib_dir)
    if environment == 'DISCOVERY':
        vib_analysis_job = job_manager.SlurmJob()
        cmd = "sbatch run.sh"
    elif environment == 'SINGLE_NODE':
        vib_analysis_job = job_manager.DefaultJob()
        cmd = "/bin/bash run.sh"
    elif environment == 'THETA':
        raise NotImplementedError
    else:
        raise NotImplementedError

    vib_analysis_job.submit(cmd)
    os.chdir(cur_dir)
