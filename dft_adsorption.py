# Sevy Harris
# 2021-11-26
# This is the main file to control all jobs


import os
import sys
import logging
import numpy as np
import yaml
from ase.io.espresso import read_espresso_out

import job_manager

START_DIR = os.path.abspath(os.path.dirname(__file__))

#################################################################
# Parse and verify input file
input_file = sys.argv[1]
with open(input_file) as f:
    input_settings = yaml.load(f, Loader=yaml.FullLoader)


BASE_DIR = os.path.abspath(os.path.dirname(input_file))
input_settings['BASE_DIR'] = BASE_DIR


if not os.path.exists(input_settings['PSEUDOS_DIR']):
    raise OSError('Must specify valid path to directory with pseudopotential files')

if not os.path.exists(input_settings['adsorbate_file']):
    raise OSError('Must specify valid path to adsorbate file')

os.makedirs(BASE_DIR, exist_ok=True)
try:
    force_recalc = input_settings['force_recalc']
except KeyError:
    force_recalc = False
force_recalc = True

try:
    job_manager_type = input_settings['job_manager']
except KeyError:
    job_manager_type = 'default'


############################################################################
# Set up logging
logfile = os.path.join(BASE_DIR, 'DFT_ADSORPTION.log')
if os.path.exists(logfile):
    os.rename(logfile, os.path.join(BASE_DIR, 'DFT_ADSORPTION.log.old'))

logfile_handle = logging.FileHandler(logfile)
stdout_handle = logging.StreamHandler(sys.stdout)
logging.basicConfig(
    format='%(asctime)s\t%(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[logfile_handle, stdout_handle]
)
logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
logger.info('Starting DFT Adsorption Calculation')


############################################################################
# Save settings to the BASE_DIR
BULK_DIR = os.path.join(BASE_DIR, 's1_bulk')
ADS_DIR = os.path.join(BASE_DIR, 's2_ads')
SLAB_DIR = os.path.join(BASE_DIR, 's3_slab')
PLACE_ADS_DIR = os.path.join(BASE_DIR, 's4_place_ads')
FINAL_SYSTEM_DIR = os.path.join(BASE_DIR, 's5_system')
input_settings['BULK_DIR'] = BULK_DIR
input_settings['ADS_DIR'] = ADS_DIR
input_settings['SLAB_DIR'] = SLAB_DIR
input_settings['PLACE_ADS_DIR'] = PLACE_ADS_DIR
input_settings['FINAL_SYSTEM_DIR'] = FINAL_SYSTEM_DIR

BULK_SLURM_JOB_PATH = os.path.abspath(os.path.join(BULK_DIR, 'run_qe.sh'))
BULK_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'bulk.py'))
BULK_PY_PATH_DST = os.path.abspath(os.path.join(BULK_DIR, 'bulk.py'))
BULK_PWO_SRC = os.path.abspath(os.path.join(BULK_DIR, 'espresso.pwo'))
BULK_PWO_DST = os.path.abspath(os.path.join(BULK_DIR, 'bulk.pwo'))
input_settings['BULK_FILE'] = BULK_PWO_DST

ADS_SLURM_JOB_PATH = os.path.abspath(os.path.join(ADS_DIR, 'run_qe.sh'))
ADS_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'adsorbate.py'))
ADS_PY_PATH_DST = os.path.abspath(os.path.join(ADS_DIR, 'adsorbate.py'))
ADS_PWO_SRC = os.path.abspath(os.path.join(ADS_DIR, 'espresso.pwo'))
ADS_PWO_DST = os.path.abspath(os.path.join(ADS_DIR, 'ads.pwo'))
input_settings['ADS_FILE'] = ADS_PWO_DST

SLAB_SLURM_JOB_PATH = os.path.abspath(os.path.join(SLAB_DIR, 'run_qe.sh'))
SLAB_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'slab.py'))
SLAB_PY_PATH_DST = os.path.abspath(os.path.join(SLAB_DIR, 'slab.py'))
SLAB_PWO_SRC = os.path.abspath(os.path.join(SLAB_DIR, 'espresso.pwo'))

# TODO check for existing slab to skip many steps
# try:
#     if os.path.exists(input_settings['SLAB_FILE']):
# except KeyError or OSError:
# SLAB_PWO_DST = os.path.abspath(os.path.join(SLAB_DIR, 'slab.pwo'))
SLAB_PWO_DST = os.path.abspath(os.path.join(SLAB_DIR, 'slab.pwo'))
input_settings['SLAB_FILE'] = SLAB_PWO_DST

settings_file = os.path.join(BASE_DIR, 'settings.yaml')
logger.info(f'Creating setting file {settings_file}')

# TODO mpi settings
logger.info(f'Settings:\n{input_settings}')
with open(settings_file, 'w') as f:
    yaml.dump(input_settings, f, sort_keys=True)

############################################################################
# Skip step 1 if there's a slab file to use
use_existing_slab = False
if os.path.exists(SLAB_PWO_DST) and not force_recalc:
    # try to read in the adsorbate
    with open(SLAB_PWO_DST, 'r') as f:
        slab_traj = list(read_espresso_out(f, index=slice(None)))
        f.seek(0)
        for line in f:
            if 'JOB DONE' in line:
                use_existing_slab = True

############################################################################
# Step 1. Compute Bulk Lattice Constant
if not use_existing_slab:
    use_existing_bulk = False
    if os.path.exists(BULK_PWO_DST) and not force_recalc:
        # try to read in the bulk lattice
        with open(BULK_PWO_DST, 'r') as f:
            bulk_traj = list(read_espresso_out(f, index=slice(None)))
            f.seek(0)
            for line in f:
                if 'JOB DONE' in line:
                    use_existing_bulk = True
                    break

    if use_existing_bulk:
        initial_energy = bulk_traj[-1].get_potential_energy()
        lattice_constant = bulk_traj[-1].cell[0][0]

        logger.info('1. Use existing bulk lattice constant')
        logger.info(f"Initial Energy: {initial_energy}")
        logger.info(f"Lattice constant: {lattice_constant}")
    else:
        logger.info('1. Compute bulk lattice constant')
        os.makedirs(BULK_DIR, exist_ok=True)

        bulk_job_script = job_manager.SlurmJobFile(full_path=BULK_SLURM_JOB_PATH)
        bulk_job_script.settings = {
            '--job-name': 'S1_BULK',
            '--partition': 'west,short',
            '--mem': '4Gb',     # max memory was 442252K
            '--time': '24:00:00',
        }
        bulk_job_script.content.append(f'cp {BULK_PY_PATH_SRC} {BULK_PY_PATH_DST}\n')
        bulk_job_script.content.append(f'python {BULK_PY_PATH_DST} {settings_file}\n')
        # bulk_job_script.content.append(f'cp {BULK_PWO_SRC} {BULK_PWO_DST}\n')
        bulk_job_script.write_file()

        # start the slurm job
        if job_manager_type == 'SLURM':
            bulk_job = job_manager.SlurmJob()
            bulk_cmd = f"sbatch {BULK_SLURM_JOB_PATH}"
        elif job_manager_type.lower() == 'default':
            bulk_job = job_manager.DefaultJob()
            bulk_cmd = f"/bin/bash {BULK_SLURM_JOB_PATH}"
        os.chdir(BULK_DIR)
        bulk_job.submit(bulk_cmd)
        os.chdir(START_DIR)


############################################################################
# Step 2. Compute Adsorbate Geometry
use_existing_ads = False
if os.path.exists(ADS_PWO_DST) and not force_recalc:
    # try to read in the adsorbate
    with open(ADS_PWO_DST, 'r') as f:
        ads_traj = list(read_espresso_out(f, index=slice(None)))
        f.seek(0)
        for line in f:
            if 'JOB DONE' in line:
                use_existing_ads = True

if use_existing_ads:
    ads_energy = ads_traj[-1].get_potential_energy()
    logger.info('2. Use existing adsorbate geometry')
    logger.info(f"Energy: {ads_energy}")

else:
    logger.info('2. Compute Adsorbate Geometry')
    os.makedirs(ADS_DIR, exist_ok=True)

    ads_job_script = job_manager.SlurmJobFile(full_path=ADS_SLURM_JOB_PATH)
    ads_job_script.settings = {
        '--job-name': 'S2_ADSORBATE',
        '--partition': 'west,short',
        '--mem': '20Gb',     # max memory TBD
        '--time': '24:00:00',
    }
    ads_job_script.content.append(f'cp {ADS_PY_PATH_SRC} {ADS_PY_PATH_DST}\n')
    ads_job_script.content.append(f'python {ADS_PY_PATH_DST} {settings_file}\n')
    # this will be handled in the python script
    # ads_job_script.content.append(f'cp {ADS_PWO_SRC} {ADS_PWO_DST}\n')
    ads_job_script.write_file()

    # start the slurm job
    if job_manager_type == 'SLURM':
        ads_job = job_manager.SlurmJob()
        ads_cmd = f"sbatch {ADS_SLURM_JOB_PATH}"
    elif job_manager_type.lower() == 'default':
        ads_job = job_manager.DefaultJob()
        ads_cmd = f"/bin/bash {ADS_SLURM_JOB_PATH}"
    os.chdir(ADS_DIR)
    ads_job.submit(ads_cmd)
    os.chdir(START_DIR)


############################################################################
# Step 3. Compute Slab Geometry
# Depends on 1
logger.info('3. Compute Slab Geometry')
os.makedirs(os.path.join(BASE_DIR, 's3_slab'), exist_ok=True)

if use_existing_slab:
    slab_energy = slab_traj[-1].get_potential_energy()
    logger.info('3. Use existing slab geometry')
    logger.info(f"Energy: {slab_energy}")

else:
    if not use_existing_bulk:
        bulk_job.wait()

    logger.info('3. Compute Slab Geometry')
    os.makedirs(SLAB_DIR, exist_ok=True)

    slab_job_script = job_manager.SlurmJobFile(full_path=SLAB_SLURM_JOB_PATH)
    slab_job_script.settings = {
        '--job-name': 'S3_SLAB',
        '--partition': 'west,short',
        '--mem': '20Gb',     # max memory TBD
        '--time': '24:00:00',
    }
    slab_job_script.content.append(f'cp {SLAB_PY_PATH_SRC} {SLAB_PY_PATH_DST}\n')
    slab_job_script.content.append(f'python {SLAB_PY_PATH_DST} {settings_file}\n')
    # slab_job_script.content.append(f'cp {SLAB_PWO_SRC} {SLAB_PWO_DST}\n')
    slab_job_script.write_file()

    # start the slurm job
    if job_manager_type == 'SLURM':
        slab_job = job_manager.SlurmJob()
        slab_cmd = f"sbatch {SLAB_SLURM_JOB_PATH}"
    elif job_manager_type.lower() == 'default':
        slab_job = job_manager.DefaultJob()
        slab_cmd = f"/bin/bash {SLAB_SLURM_JOB_PATH}"
    os.chdir(SLAB_DIR)
    slab_job.submit(slab_cmd)
    os.chdir(START_DIR)


############################################################################
# Step 4. Place Adsorbate
# Depends on 2 and 3
logger.info('4. Place Adsorbate')
os.makedirs(PLACE_ADS_DIR, exist_ok=True)

# TODO check if we need to wait on any processes
heights = np.linspace(
    input_settings['min_height_A'],
    input_settings['max_height_A'],
    input_settings['N_heights']
)
for i, height in enumerate(heights):
    os.makedirs(os.path.join(PLACE_ADS_DIR, f'run_{i}'))
    # place_ads_job_script = job_manager.SlurmJobFile(full_path=SLAB_SLURM_JOB_PATH)
    # slab_job_script.settings = {
    #     '--job-name': 'S3_SLAB',
    #     '--partition': 'west,short',
    #     '--mem': '20Gb',     # max memory TBD
    #     '--time': '24:00:00',
    # }
    # slab_job_script.content.append(f'cp {SLAB_PY_PATH_SRC} {SLAB_PY_PATH_DST}\n')
    # slab_job_script.content.append(f'python {SLAB_PY_PATH_DST} {settings_file}\n')
    # # slab_job_script.content.append(f'cp {SLAB_PWO_SRC} {SLAB_PWO_DST}\n')
    # slab_job_script.write_file()

############################################################################
# Step 5. Compute System Geometry
# Depends on 4
logger.info('5. Compute System Geometry')
# os.makedirs(os.path.join(BASE_DIR, 's5_system'), exist_ok=True)
