# Sevy Harris
# 2021-11-26
# This is the main file to control all jobs


import os
import sys
import logging
import yaml
from ase.io.espresso import read_espresso_out

import job_manager


# input variables - ask Lance for access to dftinputgen

# surface element
# adsorbate name
# surface site
# DFT functional
# energy cutoffs
# k pts sampling
# vacuum between unit cells
# force convergence

# TODO check for restart


pseudos_dir = '/home/harris.se/espresso/pseudos'  # TODO handle this with environment variable
BASE_DIR = '/scratch/harris.se/espresso/copper111/co2/'
ADS_FILE = os.path.abspath(os.path.join(os.path.dirname(__file__), 'adsorbates', 'CO2.xyz'))
BULK_DIR = os.path.join(BASE_DIR, 's1_bulk')
ADS_DIR = os.path.join(BASE_DIR, 's2_ads')
os.makedirs(BASE_DIR, exist_ok=True)
force_recalc = True
# force_recalc = False


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
# Save settings
settings_file = os.path.join(BASE_DIR, 'dft_adsorption_settings.yaml')
logger.info(f'Creating setting file {settings_file}')
settings = {
    'PSEUDOS_DIR': pseudos_dir,
    'BASE_DIR': BASE_DIR,
    'BULK_DIR': BULK_DIR,
    'ADS_DIR': ADS_DIR,
    'metal': 'Cu',
    'crystal_structure': 'fcc',
    'lattice_constant_guess': 3.6,
    'adsorbate': 'CO2',
    'adsorbate_file': ADS_FILE,
    'site': 'top',
    'dft_functional': 'default',
    'pseudopotentials': {
        'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF',
        'C': 'C.pbe-n-kjpaw_psl.1.0.0.UPF',
        'O': 'O.pbe-n-kjpaw_psl.1.0.0.UPF', 
    },
    'ecutwfc': 50,
    'ecutrho': 500,
    'kpts_bulk': (4, 4, 4),
    'kpts_slab': (4, 4, 1),
    'kpts_ads': None,
    'forc_conv_thr_eVA': 0.01,  # units of eV/A
    'forc_conv_N': 0,
    'vacuum_ads': 7.5,
}
# TODO mpi settings
logger.info(f'Settings:\n{settings}')
with open(settings_file, 'w') as f:
    yaml.dump(settings, f)


############################################################################
# Step 1. Compute Bulk Lattice Constant
bulk_file = os.path.join(BULK_DIR, 'espresso.pwo')
if os.path.exists(bulk_file) and not force_recalc:
    # try to read in the bulk lattice
    with open(bulk_file, 'r') as f:
        traj = list(read_espresso_out(f, index=slice(None)))

    initial_energy = traj[-1].get_potential_energy()
    lattice_constant = traj[-1].cell[0][0]

    logger.info('1. Use existing bulk lattice constant')
    logger.info(f"Initial Energy: {initial_energy}")
    logger.info(f"Lattice constant: {lattice_constant}")

else:
    logger.info('1. Compute bulk lattice constant')
    os.makedirs(BULK_DIR, exist_ok=True)

    # Make the slurm job
    BULK_SLURM_JOB_PATH = os.path.join(BULK_DIR, 'run_qe.sh')
    BULK_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'bulk.py'))
    BULK_PY_PATH_DST = os.path.join(BULK_DIR, 'bulk.py')
    bulk_job_script = job_manager.SlurmJobFile(full_path=BULK_SLURM_JOB_PATH)
    bulk_job_script.settings = {
        '--job-name': 'S1_BULK',
        '--partition': 'west,short',
        '--mem': '4Gb',     # max memory was 442252K
        '--time': '24:00:00',
    }
    bulk_job_script.content.append(f'cp {BULK_PY_PATH_SRC} {BULK_PY_PATH_DST}\n')
    bulk_job_script.content.append(f'python {BULK_PY_PATH_DST} {settings_file}\n')
    bulk_job_script.write_file()

    # start the slurm job
    bulk_job = job_manager.SlurmJob()
    bulk_cmd = f"sbatch {BULK_SLURM_JOB_PATH}"
    bulk_job.submit(bulk_cmd)

############################################################################
# Step 2. Compute Adsorbate Geometry
adsorbate_file = os.path.join(ADS_DIR, 'espresso.pwo')
if os.path.exists(adsorbate_file) and not force_recalc:
    # try to read in the adsorbate
    with open(adsorbate_file, 'r') as f:
        ads_traj = list(read_espresso_out(f, index=slice(None)))

    ads_energy = ads_traj[-1].get_potential_energy()
    logger.info('2. Use existing adsorbate geometry')
    logger.info(f"Energy: {ads_energy}")

else:
    logger.info('2. Compute Adsorbate Geometry')
    os.makedirs(ADS_DIR, exist_ok=True)

    # Make the slurm job
    ADS_SLURM_JOB_PATH = os.path.join(ADS_DIR, 'run_qe.sh')
    ADS_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'adsorbate.py'))
    # print(f'ADS_PY_PATH_SRC: {ADS_PY_PATH_SRC}')
    ADS_PY_PATH_DST = os.path.join(ADS_DIR, 'adsorbate.py')
    ads_job_script = job_manager.SlurmJobFile(full_path=ADS_SLURM_JOB_PATH)
    ads_job_script.settings = {
        '--job-name': 'S2_ADSORBATE',
        '--partition': 'west,short',
        '--mem': '20Gb',     # max memory TBD
        '--time': '24:00:00',
    }
    ads_job_script.content.append(f'cp {ADS_PY_PATH_SRC} {ADS_PY_PATH_DST}\n')
    ads_job_script.content.append(f'python {ADS_PY_PATH_DST} {settings_file}\n')
    ads_job_script.write_file()

    # start the slurm job
    ads_job = job_manager.SlurmJob()
    ads_cmd = f"sbatch {ADS_SLURM_JOB_PATH}"
    ads_job.submit(ads_cmd)

############################################################################
# Step 3. Compute Slab Geometry
# Depends on 1
logger.info('3. Compute Slab Geometry')
os.makedirs(os.path.join(BASE_DIR, 's3_slab'), exist_ok=True)


############################################################################
# Step 4. Place Adsorbate
# Depends on 2 and 3
logger.info('4. Place Adsorbate')
os.makedirs(os.path.join(BASE_DIR, 's4_ads_height'), exist_ok=True)


############################################################################
# Step 5. Compute System Geometry
# Depends on 4
logger.info('5. Compute System Geometry')
os.makedirs(os.path.join(BASE_DIR, 's5_system'), exist_ok=True)



