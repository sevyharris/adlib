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


def make_run_relax_script(calc_dir, nproc=32, job_name='relax_system'):
    bash_filename = os.path.join(calc_dir, 'run.sh')
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        f.write('#SBATCH --time=24:00:00\n')
        f.write(f'#SBATCH --job-name={job_name}\n')
        f.write('#SBATCH --mem=40Gb\n')
        f.write('#SBATCH --cpus-per-task=1\n')
        f.write(f'#SBATCH --ntasks={nproc}\n')
        f.write('#SBATCH --partition=short,west\n')
        f.write('module load gcc/10.1.0\n')
        f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
        f.write('module load scalapack/2.1.0-skylake\n\n')
        f.write(f'cd {calc_dir}\n')
        f.write(f'python relax_system.py\n')


# TODO actually implement this
def make_relax_script(calc_dir, lattice_constant, metal='Cu', ecutwfc=100, kpt=7, smear=0.1, nproc=32):
    """Function to make a python script to relax the slab-adsorption system
    """
    fmax = 0.01
    vacuum = 10.0

    python_file_lines = [
        "import os",
        "import sys",
        "from time import time",
        "from ase.calculators.espresso import Espresso",
        "from ase.build import fcc111",
        "from ase.constraints import FixAtoms",
        "from ase.io.ulm import InvalidULMFileError",
        "from ase.io.trajectory import Trajectory",
        "from ase.optimize import BFGS",
        "from ase.io.espresso import read_espresso_out",
        "",
        "",
        "start = time()",
        "logfile = 'ase.log'",
        "this_dir = os.path.dirname(os.path.abspath(__file__))",
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
        f"fmax = {fmax}",
        f"vacuum = {vacuum}",
        f"slab = fcc111('{metal}', size=(3, 3, 4), vacuum=vacuum, a={lattice_constant})",
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
        "slab.calc = espresso",
        "opt = BFGS(slab, logfile=logfile, trajectory='slab.traj')",
        "opt.run(fmax=fmax)",
        "",
        "# Read the energy back in",
        "energy = 0",
        "with open(os.path.join(this_dir, 'espresso.pwo'), 'r') as f:",
        "    traj = list(read_espresso_out(f, index=slice(None)))",
        "    energy = traj[-1].get_potential_energy()",
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
    os.makedirs(calc_dir, exist_ok=True)
    calc_filename = os.path.join(calc_dir, f'relax_slab.py')
    with open(calc_filename, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])


def run_relax_system(system_dir):
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(system_dir)
    relax_system_job = job_manager.SlurmJob()
    cmd = "sbatch run.sh"
    relax_system_job.submit(cmd)
    os.chdir(cur_dir)
