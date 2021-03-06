import adlib.adsorbate.convergence


adsorbate_dir = '/work/westgroup/harris.se/espresso/adlib/test/adsorbate/CO2'
adsorbate_pwo = '/work/westgroup/harris.se/espresso/qe_workflow/results/adsorbate/CO2/espresso.pwo'
adlib.adsorbate.convergence.setup_converge(adsorbate_dir, 'vacuum_converge', adsorbate_pwo=adsorbate_pwo)
adlib.adsorbate.convergence.run_converge(adsorbate_dir, 'vacuum_converge')
