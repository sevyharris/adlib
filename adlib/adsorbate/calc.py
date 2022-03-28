"""
Module for setting up and running vc-relax in Quantum Espresso

Expects/creates the following directory structure:

adsorbate_name  # something like CO2
____adsorbate.pwo  # copy of espresso.pwo
____relax_ads.py
____espresso.pwo
____run.sh
____convergence_check
________ecutwfc
____________run_qe_jobs.sh
____________run0
________________calc.py
________________espresso.pwo
____________run1
________________calc.py
________________espresso.pwo
____________runN
________________calc.py
________________espresso.pwo
________vacuum
____________run_qe_jobs.sh
____________run0
________________calc.py
________________espresso.pwo
____________run1
________________calc.py
________________espresso.pwo
____________runN
________________calc.py
________________espresso.pwo

"""

import os
import shutil
import adlib.env


environment = adlib.env.load_environment()


def setup_relax_adsorbate(adsorbate_dir, xyz_dir=None, nproc=48):
    ads_name = os.path.basename(adsorbate_dir)
    if xyz_dir is None:
        # TODO make this relative to the package and ship the code with some example adsorbate xyzs
        xyz_dir = '/work/westgroup/harris.se/espresso/qe_workflow/resources/adsorbates/'
    os.makedirs(adsorbate_dir, exist_ok=True)
    xyz_file = os.path.join(xyz_dir, f'{ads_name}.xyz')
    shutil.copy(xyz_file, adsorbate_dir)
    make_relax_ads_script(adsorbate_dir, nproc=nproc)
    make_run_relax_ads_script(adsorbate_dir, nproc=nproc)


def make_run_array(dest_dir, N_runs, job_name='ads_converge', nproc=16):
    # function to set up an array job
    bash_filename = os.path.join(dest_dir, 'run_qe_jobs.sh')
    run_i_dir = os.path.abspath(os.path.join(dest_dir, 'run_$RUN_i'))
    # write the array job file
    with open(bash_filename, 'w') as f:
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

        f.write('module load gcc/10.1.0\n')
        f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
        f.write('module load scalapack/2.1.0-skylake\n\n')

        f.write(f'cd {run_i_dir}\n')
        f.write(f'python relax_ads.py\n\n')


def make_run_relax_ads_script(calc_dir, nproc=16):
    """
    Make a shell script to run the python script
    """
    os.makedirs(calc_dir, exist_ok=True)
    ads_name = os.path.basename(calc_dir)
    bash_filename = os.path.join(calc_dir, 'run.sh')
    # write the array job file
    with open(bash_filename, 'w') as f:
        f.write('#!/bin/bash\n\n')
        if environment == 'DISCOVERY':
            f.write('#SBATCH --time=24:00:00\n')
            f.write(f'#SBATCH --job-name={ads_name}_relax' + '\n')
            f.write('#SBATCH --mem=40Gb\n')
            f.write('#SBATCH --cpus-per-task=1\n')
            f.write(f'#SBATCH --ntasks={nproc}' + '\n')
            f.write('#SBATCH --partition=short,west\n')
            f.write('module load gcc/10.1.0\n')
            f.write('module load openmpi/4.0.5-skylake-gcc10.1\n')
            f.write('module load scalapack/2.1.0-skylake\n\n')
        # f.write(f'cd {calc_dir}\n')
        f.write(f'python relax_ads.py\n')


def make_relax_ads_script(calc_dir, vacuum=10.0, ecutwfc=100, nproc=16):
    """
    Make a python script that uses ase to run quantum espresso
    """

    python_file_name = os.path.join(calc_dir, 'relax_ads.py')
    fmax = 0.01

    python_file_lines = [
        "import os",
        "import sys",
        "import glob",
        "import time",
        "from ase.calculators.espresso import Espresso",
        "from ase.optimize import BFGS",
        "from ase.io.trajectory import Trajectory",
        "from ase.io import read",
        "from ase.io.espresso import read_espresso_out",
        "from ase.io.ulm import InvalidULMFileError",
        "",
        "",
        "T = time.localtime()",
        "start = time.time()",
        "logfile = 'ase.log'",
        "with open(logfile, 'a') as f:",
        "    f.write(f'Start: {T[3]:02}:{T[4]:02}:{T[5]:02}\\n')",
        "",
        "",
        "# Assume there is only one .xyz file in the folder",
        "this_dir = os.path.dirname(os.path.abspath(__file__))",
        "files = glob.glob(os.path.join(this_dir, '*.xyz'))",
        "if len(files) > 1:",
        "    print('multiple xyzs')",
        "",
        "",
        f"fmax = {fmax}  # eV/A",
        f"vacuum = {vacuum}   # A",
        "adsorbate_file = files[0]",
        "adsorbate = read(adsorbate_file)",
        "adsorbate.center(vacuum=vacuum)",
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
        "        'degauss': 0.1,",
        f"        'ecutwfc': {ecutwfc},",
        "    },",
        "}",
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
        "}",
        "",
        "",
        "# Restart if available",
        "traj_file = 'ads.traj'",
        "if os.path.exists(traj_file):",
        "    try:",
        "        traj = Trajectory(traj_file)",
        "        if len(traj) > 0:",
        "            adsorbate = traj[-1]",
        "    except InvalidULMFileError:",
        "        # traj file is empty. delete it.",
        "        os.remove(traj_file)",
        "",
        f"command = f'mpirun -np {nproc} " + "{pw_executable} -in PREFIX.pwi > PREFIX.pwo'",
        "print(command)",
        "",
        "espresso = Espresso(",
        "    command=command,",
        "    pseudopotentials=pseudopotentials,",
        "    tstress=True,",
        "    tprnfor=True,",
        "    kpts=(1, 1, 1),",
        "    pseudo_dir=os.environ['PSEUDO_DIR'],",
        "    input_data=espresso_settings,",
        ")",
        "",
        "",
        "opt = BFGS(adsorbate, trajectory=traj_file, logfile=logfile)",
        "adsorbate.calc = espresso",
        "opt.run(fmax)",
        "",
        "",
        "energy = 0",
        "with open(os.path.join(this_dir, 'espresso.pwo'), 'r') as f:",
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
    with open(python_file_name, 'w') as f:
        f.writelines([line + '\n' for line in python_file_lines])


def run_relax_adsorbate(adsorbate_dir):
    import job_manager
    cur_dir = os.getcwd()
    os.chdir(adsorbate_dir)
    if environment == 'DISCOVERY':
        relax_ads_job = job_manager.SlurmJob()
        cmd = "sbatch run.sh"
    elif environment == 'SINGLE_NODE':
        relax_ads_job = job_manager.DefaultJob()
        cmd = "/bin/bash run.sh"
    elif environment == 'THETA':
        raise NotImplementedError
    else:
        raise NotImplementedError
    relax_ads_job.submit(cmd)
    os.chdir(cur_dir)
