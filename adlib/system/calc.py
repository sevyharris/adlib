"""
Module for setting up and running system energy calculations in Quantum Espresso

system
____system.pwo
____relax_system.py
____espresso.pwo
____convergence_check
________kpts
____________run_qe_jobs.sh
____________run0000
________________calc.py
________________espresso.pwo
____________run0001
________________calc.py
________________espresso.pwo
____________run000N
________________calc.py
________________espresso.pwo
________ecutwfc
____________run_qe_jobs.sh
____________run0000
________________calc.py
________________espresso.pwo
____________run0001
________________calc.py
________________espresso.pwo
____________run000N
________________calc.py
________________espresso.pwo
________smear
____________run_qe_jobs.sh
____________run0000
________________calc.py
________________espresso.pwo
____________run0001
________________calc.py
________________espresso.pwo
____________run000N
________________calc.py
________________espresso.pwo
"""

import os
import adlib.env


environment = adlib.env.load_environment()


def make_run_relax_script(calc_dir, nproc=48, job_name='relax_system'):
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
            f.write('#SBATCH --partition=short\n')
            f.write('#SBATCH --constraint=cascadelake\n\n')
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')
        # f.write(f'cd {calc_dir}\n')
        f.write(f'python relax_system.py\n')


def make_relax_script(calc_dir, ecutwfc=50, kpt=5, smear=0.1, nproc=48, low_mixing_beta=False):
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
        "import numpy as np",
        "from ase.build import bulk, fcc111, add_adsorbate",
        "from ase.io import read, write",
        "from ase import Atoms",
        "from ase.optimize import BFGS",
        "from ase.calculators.espresso import Espresso, EspressoProfile",
        "from ase.constraints import FixAtoms",
        "from ase.io.espresso import read_espresso_out",
        "from ase.io.trajectory import Trajectory",
        "from ase.io.ulm import InvalidULMFileError",
        "import time",
        "",
        "",
        "def place_adsorbate_top(metal_slab, adsorbate, height=1.0, ads_index=0):",
        "    # remove the adsorbate cell",
        "    adsorbate.cell = [0, 0, 0]",
        "    top_layer_z = np.max([pos[2] for pos in metal_slab.get_positions()])",
        "    for i, pos in enumerate(metal_slab.get_positions()):",
        "        if pos[2] == top_layer_z:",
        "            # place the atom height above here",
        "            ads_origin = pos + [0, 0, height]",
        "            pos_difference = ads_origin - adsorbate.get_positions()[ads_index]",
        "            adsorbate.translate(pos_difference)",
        "            metal_slab += adsorbate",
        "            print(f'Using {adsorbate[ads_index]} for binding atom')",
        "            return",
        "    print('Failed to place adsorbate')",
        "    exit(-1)",
        "",
        "",
        "# TODO accept height guess",
        "# if len(sys.argv) < 2:",
        "#     raise IndexError('Must specify starting height for fine system run')",
        "# height = float(sys.argv[1])",
        "height = 3.2",
        f"fmax = {fmax}",
        "T = time.localtime()",
        "start = time.time()",
        "",
        "logfile = 'ase.log'",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Start: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "",
        "slab_file = 'slab.pwo'",
        "with open(slab_file, 'r') as f:",
        "    traj = list(read_espresso_out(f, index=slice(None)))",
        "metal_slab = traj[-1]",
        "",
        "adsorbate_file = 'adsorbate.pwo'",
        "with open(adsorbate_file, 'r') as f:",
        "    traj = list(read_espresso_out(f, index=slice(None)))",
        "adsorbate = traj[-1]",
        "",
        "# Fix the bottom two layers",
        "bottom_layer = []",
        "second_layer = []",
        "fixed_indices = []",
        "# Round to the nearest 0.1 Angstrom in determining layers",
        "z_values = list(set([np.round(pos[2], 1) for pos in metal_slab.get_positions()]))",
        "z_values.sort()",
        "",
        "",
        "for i, pos in enumerate(metal_slab.get_positions()):",
        "    if np.round(pos[2], 1) == z_values[0]:",
        "        bottom_layer.append(metal_slab[i].index)",
        "    if np.round(pos[2], 1) == z_values[1]:",
        "        second_layer.append(metal_slab[i].index)",
        "fixed_indicies = bottom_layer + second_layer",
        "fix_bottom_layers = FixAtoms(indices=fixed_indicies)",
        "metal_slab.set_constraint(fix_bottom_layers)",
        "",
        "",
        "# place the adsorbate",
        "element_priority = ['C', 'O', 'H']  # where does N fit?",
        "bond_atom_index = -1",
        "for element in element_priority:",
        "    for atom in adsorbate:",
        "        if atom.symbol == element:",
        "            bond_atom_index = atom.index",
        "            break",
        "    if bond_atom_index > -1:",
        "        break",
        "print(f'bond atom is {adsorbate[bond_atom_index]}')",
        "",
        "# 'ontop', 'bridge', 'fcc', 'hcp'",
        "# https://wiki.fysik.dtu.dk/ase/ase/build/surface.html#ase.build.add_adsorbate",
        "# add_adsorbate(metal_slab, adsorbate, height=height, position='ontop', mol_index=bond_atom_index)",
        "# can't use the add_adsorbate function if you load the atoms from a pwo file...",
        "place_adsorbate_top(metal_slab, adsorbate, height=height, ads_index=bond_atom_index)",
        "",
        "",
        "# Restart if available",
        "traj_file = 'system.traj'",
        "restart = False",
        "if os.path.exists(traj_file):",
        "    try:",
        "        traj = Trajectory(traj_file)",
        "        if len(traj) > 0:",
        "            metal_slab = traj[-1]",
        "            restart = True",
        "            with open(logfile, 'a') as f:",
        "                f.write(f'Restarting from trajectory file {traj_file}\\n')",
        "    except InvalidULMFileError:",
        "        # traj file is empty. delete it.",
        "        os.remove(traj_file)",
        "",
        "if not restart:",
        "    write(f'initial_system.xyz', metal_slab)",
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
        "}",
        "",
        "pw_executable = os.environ['PW_EXECUTABLE']",
        f"command = f'mpirun -np {nproc} " + "{pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "profile = EspressoProfile(argv=command.split())",
        "calc = Espresso(",
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
        "metal_slab.calc = calc",
        "opt = BFGS(metal_slab, logfile=logfile, trajectory='system.traj')",
        "opt.run(fmax=fmax)",
        "",
        "",
        "# Read the energy back in",
        "energy = 0",
        "with open('espresso.pwo', 'r') as f:",
        "    traj = list(read_espresso_out(f, index=slice(None)))",
        "    energy = traj[-1].get_potential_energy()",
        "",
        "",
        "T = time.localtime()",
        "end = time.time()",
        "duration = end - start",
        "",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Energy: {energy}\\n')",
        "    f.write(f'End: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
    ]

    os.makedirs(calc_dir, exist_ok=True)
    calc_filename = os.path.join(calc_dir, f'relax_system.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])


def run_relax_system(system_dir):
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(system_dir)
    if environment == 'DISCOVERY':
        relax_system_job = job_manager.SlurmJob()
        cmd = "sbatch run.sh"
    elif environment == 'SINGLE_NODE':
        relax_system_job = job_manager.DefaultJob()
        cmd = "/bin/bash run.sh"
    elif environment == 'THETA':
        raise NotImplementedError
    else:
        raise NotImplementedError

    relax_system_job.submit(cmd)
    os.chdir(cur_dir)


def make_scf_script(calc_dir, ecutwfc=60, kpt=5, smear=0.1, nproc=48, low_mixing_beta=False):
    """Function to make a python script to run scf given an input geometry
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
        "import numpy as np",
        "from ase.build import bulk, fcc111, add_adsorbate",
        "from ase.io import read, write",
        "from ase import Atoms",
        "from ase.optimize import BFGS",
        "from ase.calculators.espresso import Espresso, EspressoProfile",
        "from ase.constraints import FixAtoms",
        "from ase.io.espresso import read_espresso_out",
        "from ase.io.trajectory import Trajectory",
        "from ase.io.ulm import InvalidULMFileError",
        "import time",
        "",
        "",
        "",
        "T = time.localtime()",
        "start = time.time()",
        "",
        "logfile = 'ase.log'",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Start: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "",
        "input_geometry_file = 'input_geometry.pwo'",
        "with open(input_geometry_file, 'r') as f:",
        "    traj = list(read_espresso_out(f, index=slice(None)))",
        "system = traj[-1]",
        "",
        "",
        "espresso_settings = {",
        "    'control': {",
        "        'verbosity': 'high',",
        "        'calculation': 'scf',",
        "        'disk_io': 'none',",
        "        'max_seconds': 84600,  # 23.5 hours",
        "        'restart_mode': 'restart',",
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
        "}",
        "",
        "pw_executable = os.environ['PW_EXECUTABLE']",
        f"command = f'mpirun -np {nproc} " + "{pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "profile = EspressoProfile(argv=command.split())",
        "calc = Espresso(",
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
        "system.calc = calc",
        "energy = system.get_potential_energy()",
        "",
        "",
        "T = time.localtime()",
        "end = time.time()",
        "duration = end - start",
        "",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Energy: {energy}\\n')",
        "    f.write(f'End: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
    ]

    os.makedirs(calc_dir, exist_ok=True)
    calc_filename = os.path.join(calc_dir, f'run_scf.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])

def make_run_scf_script(calc_dir, nproc=48, job_name='run_scf'):
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
            f.write('#SBATCH --partition=short\n')
            f.write('#SBATCH --constraint=cascadelake\n\n')
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')
        # f.write(f'cd {calc_dir}\n')
        f.write(f'python run_scf.py\n')


def run_scf(calc_dir):
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(calc_dir)
    if environment == 'DISCOVERY':
        scf_job = job_manager.SlurmJob()
        cmd = "sbatch run.sh"
    elif environment == 'SINGLE_NODE':
        scf_job = job_manager.DefaultJob()
        cmd = "/bin/bash run.sh"
    elif environment == 'THETA':
        raise NotImplementedError
    else:
        raise NotImplementedError

    scf_job.submit(cmd)
    os.chdir(cur_dir)
