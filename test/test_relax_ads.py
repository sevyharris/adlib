import adlib.adsorbate.calc


adsorbate_dir = '/work/westgroup/harris.se/espresso/adlib/test/adsorbate/CO2'

adlib.adsorbate.calc.setup_relax_adsorbate(adsorbate_dir)
adlib.adsorbate.calc.run_relax_adsorbate(adsorbate_dir)
