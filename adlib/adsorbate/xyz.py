"""
Module for making initial guesses of adsorbate geometry and saving as xyz
"""
import os

from ase.io import write
from ase import Atoms
from ase.collections import g2


def make_xyz(ads_name, base_dir='.'):
    try:
        ads = g2[ads_name]
        ads.rotate(90, 'y')  # we want the adsorbate to be flat in the xy plane
        write(os.path.join(base_dir, f'{ads_name}.xyz'), ads)
    except KeyError:
        raise NotImplementedError
