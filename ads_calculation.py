# Sevy Harris
# 2021-12-27
# Module to set up jobs for dft binding energy calculation


import os
import sys
import logging
import numpy as np
import yaml
from ase.io.espresso import read_espresso_out

import job_manager


class AdsorptionCalculation():
    def __init__(self, input_file):
        self.START_DIR = os.path.abspath(os.path.dirname(__file__))
        ############################################################################
        # Parse and verify input file
        self.input_file = input_file
        with open(self.input_file) as f:
            self.input_settings = yaml.load(f, Loader=yaml.FullLoader)

        self.BASE_DIR = os.path.abspath(os.path.dirname(self.input_file))
        self.input_settings['BASE_DIR'] = self.BASE_DIR

        if not os.path.exists(self.input_settings['PSEUDOS_DIR']):
            raise OSError('Must specify valid path to directory with pseudopotential files')

        if not os.path.exists(self.input_settings['adsorbate_file']):
            raise OSError('Must specify valid path to adsorbate file')

        os.makedirs(self.BASE_DIR, exist_ok=True)

        try:
            self.job_manager_type = self.input_settings['job_manager']
        except KeyError:
            self.job_manager_type = 'default'

        ############################################################################
        # Set up logging
        self.logfile = os.path.join(self.BASE_DIR, 'DFT_ADSORPTION.log')
        if os.path.exists(self.logfile):
            os.rename(self.logfile, os.path.join(self.BASE_DIR, 'DFT_ADSORPTION.log.old'))

        self.logfile_handle = logging.FileHandler(self.logfile)
        stdout_handle = logging.StreamHandler(sys.stdout)
        logging.basicConfig(
            format='%(asctime)s\t%(message)s',
            level=logging.INFO,
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[self.logfile_handle, stdout_handle]
        )
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)
        self.logger.info('Starting DFT Adsorption Calculation')

        ############################################################################
        # Save settings to the BASE_DIR
        self.BULK_DIR = os.path.join(self.BASE_DIR, 's1_bulk')
        self.ADS_DIR = os.path.join(self.BASE_DIR, 's2_ads')
        self.SLAB_DIR = os.path.join(self.BASE_DIR, 's3_slab')
        self.PLACE_ADS_DIR = os.path.join(self.BASE_DIR, 's4_place_ads')
        self.FINAL_SYSTEM_DIR = os.path.join(self.BASE_DIR, 's5_system')
        self.input_settings['BULK_DIR'] = self.BULK_DIR
        self.input_settings['ADS_DIR'] = self.ADS_DIR
        self.input_settings['SLAB_DIR'] = self.SLAB_DIR
        self.input_settings['PLACE_ADS_DIR'] = self.PLACE_ADS_DIR
        self.input_settings['FINAL_SYSTEM_DIR'] = self.FINAL_SYSTEM_DIR

        self.BULK_SLURM_JOB_PATH = os.path.abspath(os.path.join(self.BULK_DIR, 'run_qe.sh'))
        self.BULK_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'bulk.py'))
        self.BULK_PY_PATH_DST = os.path.abspath(os.path.join(self.BULK_DIR, 'bulk.py'))
        self.BULK_PWO_SRC = os.path.abspath(os.path.join(self.BULK_DIR, 'espresso.pwo'))
        self.BULK_PWO_DST = os.path.abspath(os.path.join(self.BULK_DIR, 'bulk.pwo'))
        self.input_settings['BULK_FILE'] = self.BULK_PWO_DST

        self.ADS_SLURM_JOB_PATH = os.path.abspath(os.path.join(self.ADS_DIR, 'run_qe.sh'))
        self.ADS_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'adsorbate.py'))
        self.ADS_PY_PATH_DST = os.path.abspath(os.path.join(self.ADS_DIR, 'adsorbate.py'))
        self.ADS_PWO_SRC = os.path.abspath(os.path.join(self.ADS_DIR, 'espresso.pwo'))
        self.ADS_PWO_DST = os.path.abspath(os.path.join(self.ADS_DIR, 'ads.pwo'))
        self.input_settings['ADS_FILE'] = self.ADS_PWO_DST

        self.SLAB_SLURM_JOB_PATH = os.path.abspath(os.path.join(self.SLAB_DIR, 'run_qe.sh'))
        self.SLAB_PY_PATH_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), 'basepys', 'slab.py'))
        self.SLAB_PY_PATH_DST = os.path.abspath(os.path.join(self.SLAB_DIR, 'slab.py'))
        self.SLAB_PWO_SRC = os.path.abspath(os.path.join(self.SLAB_DIR, 'espresso.pwo'))

        # TODO check for existing slab to skip many steps
        # try:
        #     if os.path.exists(input_settings['SLAB_FILE']):
        # except KeyError or OSError:
        # SLAB_PWO_DST = os.path.abspath(os.path.join(SLAB_DIR, 'slab.pwo'))
        self.SLAB_PWO_DST = os.path.abspath(os.path.join(self.SLAB_DIR, 'slab.pwo'))
        self.input_settings['SLAB_FILE'] = self.SLAB_PWO_DST

        self.settings_file = os.path.join(self.BASE_DIR, 'settings.yaml')
        self.logger.info(f'Creating setting file {self.settings_file}')

        # TODO mpi settings
        self.logger.info(f'Settings:\n{self.input_settings}')
        with open(self.settings_file, 'w') as f:
            yaml.dump(self.input_settings, f, sort_keys=True)

    def calc_bulk(self, force_recalc=False):
        """
        User's responsibility to check for existing slab
        """
        # Step 1. Compute Bulk Lattice Constant
        use_existing_bulk = False
        if os.path.exists(self.BULK_PWO_DST) and not force_recalc:
            # try to read in the bulk lattice
            with open(self.BULK_PWO_DST, 'r') as f:
                bulk_traj = list(read_espresso_out(f, index=slice(None)))
                f.seek(0)
                for line in f:
                    if 'JOB DONE' in line:
                        use_existing_bulk = True
                        break

        if use_existing_bulk:
            initial_energy = bulk_traj[-1].get_potential_energy()
            lattice_constant = bulk_traj[-1].cell[0][0]

            self.logger.info('1. Use existing bulk lattice constant')
            self.logger.info(f"Initial Energy: {initial_energy}")
            self.logger.info(f"Lattice constant: {lattice_constant}")
        else:
            self.logger.info('1. Compute bulk lattice constant')
            os.makedirs(self.BULK_DIR, exist_ok=True)

            bulk_job_script = job_manager.SlurmJobFile(full_path=self.BULK_SLURM_JOB_PATH)
            bulk_job_script.settings = {
                '--job-name': 'S1_BULK',
                '--partition': 'west,short',
                '--mem': '4Gb',     # max memory was 442252K
                '--time': '24:00:00',
            }
            bulk_job_script.content.append(f'cp {self.BULK_PY_PATH_SRC} {self.BULK_PY_PATH_DST}\n')
            bulk_job_script.content.append(f'python {self.BULK_PY_PATH_DST} {self.settings_file}\n')
            bulk_job_script.write_file()

            # start the slurm job
            if self.job_manager_type == 'SLURM':
                bulk_job = job_manager.SlurmJob()
                bulk_cmd = f"sbatch {self.BULK_SLURM_JOB_PATH}"
            elif self.job_manager_type.lower() == 'default':
                bulk_job = job_manager.DefaultJob()
                bulk_cmd = f"/bin/bash {self.BULK_SLURM_JOB_PATH}"
            os.chdir(self.BULK_DIR)
            bulk_job.submit(bulk_cmd)
            os.chdir(self.START_DIR)

    def calc_adsorbate(self, force_recalc=False):
        use_existing_ads = False
        if os.path.exists(self.ADS_PWO_DST) and not force_recalc:
            # try to read in the adsorbate
            with open(self.ADS_PWO_DST, 'r') as f:
                ads_traj = list(read_espresso_out(f, index=slice(None)))
                f.seek(0)
                for line in f:
                    if 'JOB DONE' in line:
                        use_existing_ads = True

        if use_existing_ads:
            ads_energy = ads_traj[-1].get_potential_energy()
            self.logger.info('2. Use existing adsorbate geometry')
            self.logger.info(f"Energy: {ads_energy}")

        else:
            self.logger.info('2. Compute Adsorbate Geometry')
            os.makedirs(self.ADS_DIR, exist_ok=True)

            ads_job_script = job_manager.SlurmJobFile(full_path=self.ADS_SLURM_JOB_PATH)
            ads_job_script.settings = {
                '--job-name': 'S2_ADSORBATE',
                '--partition': 'west,short',
                '--mem': '20Gb',     # max memory TBD
                '--time': '24:00:00',
            }
            ads_job_script.content.append(f'cp {self.ADS_PY_PATH_SRC} {self.ADS_PY_PATH_DST}\n')
            ads_job_script.content.append(f'python {self.ADS_PY_PATH_DST} {self.settings_file}\n')
            ads_job_script.write_file()

            # start the slurm job
            if self.job_manager_type == 'SLURM':
                ads_job = job_manager.SlurmJob()
                ads_cmd = f"sbatch {self.ADS_SLURM_JOB_PATH}"
            elif self.job_manager_type.lower() == 'default':
                ads_job = job_manager.DefaultJob()
                ads_cmd = f"/bin/bash {self.ADS_SLURM_JOB_PATH}"
            os.chdir(self.ADS_DIR)
            ads_job.submit(ads_cmd)
            os.chdir(self.START_DIR)

    def calc_slab(self, force_recalc=False):
        self.logger.info('3. Compute Slab Geometry')
        os.makedirs(self.SLAB_DIR, exist_ok=True)

        # Skip step 1 if there's a slab file to use
        use_existing_slab = False
        if os.path.exists(self.SLAB_PWO_DST) and not force_recalc:
            # try to read in the adsorbate
            with open(self.SLAB_PWO_DST, 'r') as f:
                slab_traj = list(read_espresso_out(f, index=slice(None)))
                f.seek(0)
                for line in f:
                    if 'JOB DONE' in line:
                        use_existing_slab = True

        if use_existing_slab:
            slab_energy = slab_traj[-1].get_potential_energy()
            self.logger.info('3. Use existing slab geometry')
            self.logger.info(f"Energy: {slab_energy}")

        else:
            # TODO wait on job completion
            # if not use_existing_bulk:
            #     bulk_job.wait()

            self.logger.info('3. Compute Slab Geometry')
            os.makedirs(self.SLAB_DIR, exist_ok=True)

            slab_job_script = job_manager.SlurmJobFile(full_path=self.SLAB_SLURM_JOB_PATH)
            slab_job_script.settings = {
                '--job-name': 'S3_SLAB',
                '--partition': 'west,short',
                '--mem': '20Gb',     # max memory TBD
                '--time': '24:00:00',
            }
            slab_job_script.content.append(f'cp {self.SLAB_PY_PATH_SRC} {self.SLAB_PY_PATH_DST}\n')
            slab_job_script.content.append(f'python {self.SLAB_PY_PATH_DST} {self.settings_file}\n')
            slab_job_script.write_file()

            # start the slurm job
            if self.job_manager_type == 'SLURM':
                slab_job = job_manager.SlurmJob()
                slab_cmd = f"sbatch {self.SLAB_SLURM_JOB_PATH}"
            elif self.job_manager_type.lower() == 'default':
                slab_job = job_manager.DefaultJob()
                slab_cmd = f"/bin/bash {self.SLAB_SLURM_JOB_PATH}"
            os.chdir(self.SLAB_DIR)
            slab_job.submit(slab_cmd)
            os.chdir(self.START_DIR)

    def bulk_calc_completed(self):
        pass
        # maybe try to read in the objects from the pwo files...
        # either the calculation is in progress, hasn't been done, or was done

# ############################################################################
# # Step 4. Place Adsorbate
# # Depends on 2 and 3
# logger.info('4. Place Adsorbate')
# os.makedirs(PLACE_ADS_DIR, exist_ok=True)

# # TODO check if we need to wait on any processes
# heights = np.linspace(
#     input_settings['min_height_A'],
#     input_settings['max_height_A'],
#     input_settings['N_heights']
# )
# for i, height in enumerate(heights):
#     os.makedirs(os.path.join(PLACE_ADS_DIR, f'run_{i}'))
#     # place_ads_job_script = job_manager.SlurmJobFile(full_path=SLAB_SLURM_JOB_PATH)
#     # slab_job_script.settings = {
#     #     '--job-name': 'S3_SLAB',
#     #     '--partition': 'west,short',
#     #     '--mem': '20Gb',     # max memory TBD
#     #     '--time': '24:00:00',
#     # }
#     # slab_job_script.content.append(f'cp {SLAB_PY_PATH_SRC} {SLAB_PY_PATH_DST}\n')
#     # slab_job_script.content.append(f'python {SLAB_PY_PATH_DST} {settings_file}\n')
#     # # slab_job_script.content.append(f'cp {SLAB_PWO_SRC} {SLAB_PWO_DST}\n')
#     # slab_job_script.write_file()

# ############################################################################
# # Step 5. Compute System Geometry
# # Depends on 4
# logger.info('5. Compute System Geometry')
# # os.makedirs(os.path.join(BASE_DIR, 's5_system'), exist_ok=True)
