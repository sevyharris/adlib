# Sevy Harris
# 2021-11-26
# This is the main file to control all jobs

import os
import sys
import logging
from datetime import datetime

# input variables - ask Lance for access to dftinputgen

# surface element
# adsorbate name
# surface site
# DFT functional
# energy cutoffs
# k pts sampling
# vacuum between unit cells
# force convergence


basedir = '/scratch/harris.se/espresso/copper111/co2/'
os.makedirs(basedir, exist_ok=True)
logfile = os.path.join(basedir, 'DFT_ADSORPTION.log')
if os.path.exists(logfile):
    os.rename(logfile, os.path.join(basedir, 'DFT_ADSORPTION.log.old'))

logfile_handle = logging.FileHandler(logfile)
stdout_handle = logging.StreamHandler(sys.stdout)
#logger.addHandler(logfile_handle)
#logger.addHandler(stdout_handle)

logging.basicConfig(
    format='%(asctime)s\t%(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[logfile_handle, stdout_handle]
)

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)

logger.info('Starting DFT Adsorption Calculation')


