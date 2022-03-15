import os
from shutil import copy
import adlib.system.calc

# Needs the slab.pwo and adsorbate.pwo in the directory
system_dir = '/work/westgroup/harris.se/espresso/adlib/test/system'

# replace with paths to slab and adsorbate pwo's
copy('/work/westgroup/harris.se/espresso/qe_workflow/resources/slab.pwo', system_dir)
copy('/work/westgroup/harris.se/espresso/qe_workflow/results/adsorbate/CO2/espresso.pwo', os.path.join(system_dir, 'adsorbate.pwo'))

adlib.system.calc.make_relax_script(system_dir, 3.6)
adlib.system.calc.make_run_relax_script(system_dir)
adlib.system.calc.run_relax_system(system_dir)
