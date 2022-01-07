import os
import glob
import shutil


nodes_per_python = 16
#BASE_DIR = os.getcwd()
BASE_DIR = "/home/sevy/espresso/adsorbates/"
input_ads_dir = "/home/sevy/espresso/dft_adsorption/adsorbates"
xyzs = glob.glob(os.path.join(input_ads_dir, '*.xyz'))
# or specify the xyzs manually

adsorbate_calculation_file = '/home/sevy/espresso/dft_adsorption/basepys/adsorbate.py'


python_files = []
for xyz in xyzs:
    # Make the adsorbate directory
    base_name = os.path.split(xyz)[1]
    ads_name = os.path.splitext(base_name)[0]
    # print(ads_name)
    ads_dir = os.path.join(BASE_DIR, ads_name)
    os.makedirs(ads_dir, exist_ok=True)

    
    # Copy the corresponding .xyz file into the folder
    shutil.copy(xyz, ads_dir)

    # Copy the ase calculation script into that directory
    shutil.copy(adsorbate_calculation_file, ads_dir)
    base_python = os.path.split(adsorbate_calculation_file)[1]
    python_files.append(os.path.join(ads_dir, base_python))


total_nodes = nodes_per_python * len(python_files)

with open('run.sh', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('#COBALT -q debug-cache-quad\n')
    f.write('#COBALT -A catalysis_aesp\n')
    f.write(f'#COBALT -n {total_nodes}\n')
    f.write('#COBALT -t 3:00:00\n')
    f.write('#COBALT -O pw\n\n')
    
    f.write('module load miniconda-3\n\n')
    
    for python_file in python_files:
        adsorbate_dir = os.path.dirname(python_file)
        f.write(f'cd {adsorbate_dir}\n')
        f.write(f'python {python_file} &\n')
        f.write(f'cd ..\n')
        f.write(f'sleep 5\n\n')

    f.write(f'wait\n')
