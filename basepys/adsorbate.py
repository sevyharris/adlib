# Module to compute isolated adsorbate geometry/energy
# Sevy Harris
# 2022-01-05

import glob
import os
import shutil
import sys
import socket
from shutil import copyfile
from time import time
from ase.calculators.espresso import Espresso
from ase.calculators.socketio import SocketIOCalculator
from ase import Atoms
from ase.io import read
from ase.optimize import BFGS
from ase.io.trajectory import Trajectory
from ase.io.ulm import InvalidULMFileError


start = time()
logfile = 'ase.log'


# Assume there is only one .xyz file in the folder
this_dir = os.path.dirname(os.path.abspath(__file__))

files = glob.glob(os.path.join(this_dir, "*.xyz"))
if len(files) > 1:
    print("multiple xyz's")

fmax = 0.001  # eV/A
vacuum = 15.0   # A
adsorbate_file = files[0]
adsorbate = read(adsorbate_file)
adsorbate.center(vacuum=vacuum)


espresso_settings = {
    'control': {
        'disk_io': 'none',
        'calculation': 'scf',
    },
    'system': {
        'input_dft': 'BEEF-VDW',
        'occupations': 'smearing',
        'degauss': 0.1,
        'ecutwfc': 100,
        #'ecutrho': 500,
    },
}



hostname = socket.gethostname()
# port = 31415  # the default port

min_port = 10000
max_port = 40000
seed = time() * 100000
port = int(seed % (max_port - min_port) + min_port)

pw_executable = os.environ['PW_EXECUTABLE']
pseudopotentials = {
    'C': 'C_ONCV_PBE-1.2.upf',
    'Cu': 'Cu_ONCV_PBE-1.2.upf',
    'O': 'O_ONCV_PBE-1.2.upf',
    'N': 'N_ONCV_PBE-1.2.upf',
    'H': 'H_ONCV_PBE-1.2.upf',
}


# Restart if available
traj_file = 'ads.traj'
if os.path.exists(traj_file):
    try:
        traj = Trajectory(traj_file)
        if len(traj) > 0:
            adsorbate = traj[-1]
    except InvalidULMFileError:
        # traj file is empty. delete it.
        os.remove(traj_file)


#aprun -n 4 -N 1 pw.x -nk 4 --ipi {host}:{port} --in PREFIX.pwi > PREFIX.out
command = f'aprun -n 16 {pw_executable} -in PREFIX.pwi --ipi {hostname}:{port} -nk 4 > PREFIX.pwo'
print(command)
espresso = Espresso(
    command=command,
    pseudopotentials=pseudopotentials,
    tstress=True,
    tprnfor=True,
    kpts=(4, 4, 4),
    pseudo_dir='/home/sharris1585/espresso/pseudos/',
    input_data=espresso_settings,
)

opt = BFGS(adsorbate, trajectory=traj_file, logfile=logfile)
with SocketIOCalculator(espresso, log=sys.stdout, port=port) as calc:
    adsorbate.calc = calc
    opt.run(fmax)


# copy the file
shutil.copyfile('espresso.pwo', 'adsorbate.pwo')

end = time()
duration = end - start

with open(logfile, 'a') as f:
    f.write(f'Completed in {duration} seconds\n')

