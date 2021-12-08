#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --job-name=DFT_ADSORPTION
#SBATCH --partition=west,short
#SBATCH --mem=20Gb
#SBATCH --ntasks=4

python dft_adsorption.py co2_cu_beefvdw.yaml


