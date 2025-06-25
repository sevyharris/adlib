"""
Module for setting up and running slab energy calculations in Quantum Espresso

slab
____slab.pwo
____relax_slab.py
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


def make_run_relax_script(calc_dir, nproc=48, job_name='relax_slab'):
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
        f.write(f'python relax_slab.py\n')


def make_relax_script(calc_dir, lattice_constant, metal='Cu', ecutwfc=60, kpt=5, smear=0.1, nproc=48, slab_size=(3, 3, 3)):
    """Function to make a python script to relax the slab
    """
    fmax = 0.01
    vacuum = 7.5

    python_file_lines = [
        "import os",
        "import sys",
        "import time",
        "from ase.calculators.espresso import Espresso, EspressoProfile",
        "from ase.build import fcc111",
        "from ase.constraints import FixAtoms",
        "from ase.io.ulm import InvalidULMFileError",
        "from ase.io.trajectory import Trajectory",
        "from ase.optimize import BFGS",
        "from ase.io.espresso import read_espresso_out",
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
        "    },",
        "}",
        "",
        "",
        f"fmax = {fmax}",
        f"vacuum = {vacuum}",
        f"slab = fcc111('{metal}', size={slab_size}, vacuum=vacuum, a={lattice_constant})",
        "",
        "",
        "# Fix the bottom two layers",
        "bottom_layer = []",
        "second_layer = []",
        "fixed_indices = []",
        "z_values = list(set([pos[2] for pos in slab.get_positions()]))",
        "z_values.sort()",
        "",
        "",
        "for i, pos in enumerate(slab.get_positions()):",
        "    if pos[2] == z_values[0]:",
        "       bottom_layer.append(slab[i].index)",
        "    if pos[2] == z_values[1]:",
        "       second_layer.append(slab[i].index)",
        "fixed_indicies = bottom_layer + second_layer",
        "fix_bottom_layers = FixAtoms(indices=fixed_indicies)",
        "slab.set_constraint(fix_bottom_layers)",
        "",
        "",
        "# Restart if available",
        "traj_file = 'slab.traj'",
        "if os.path.exists(traj_file):",
        "    try:",
        "        traj = Trajectory(traj_file)",
        "        if len(traj) > 0:",
        "            slab = traj[-1]",
        "    except InvalidULMFileError:",
        "        # traj file is empty. delete it.",
        "        os.remove(traj_file)",
        "",
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
        f"    kpts=({kpt}, {kpt}, 1),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "slab.calc = espresso",
        "opt = BFGS(slab, logfile=logfile, trajectory='slab.traj')",
        "opt.run(fmax=fmax)",
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
        "    f.write(f'Energy: {energy} eV\\n')",
        "    f.write(f'End: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "    f.write(f'Completed in {duration} seconds\\n')",
        "",
    ]
    os.makedirs(calc_dir, exist_ok=True)
    calc_filename = os.path.join(calc_dir, f'relax_slab.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])


def run_relax_slab(slab_dir):
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(slab_dir)
    if environment == 'DISCOVERY':
        relax_slab_job = job_manager.SlurmJob()
        cmd = "sbatch run.sh"
    elif environment == 'SINGLE_NODE':
        relax_slab_job = job_manager.DefaultJob()
        cmd = "/bin/bash run.sh"
    elif environment == 'THETA':
        raise NotImplementedError
    else:
        raise NotImplementedError

    relax_slab_job.submit(cmd)
    os.chdir(cur_dir)
