import adlib.adsorbate.xyz

xyz_dir = '/work/westgroup/harris.se/espresso/adlib/test/adsorbate/xyzs'

adsorbates = [
    'CO2',
    'CO',
    'O',
    'H',
    'OH',
    'H2O',
    'CH4',
    'C2H4',
    'N2',
    'O2',
    'H2',
]
for ads in adsorbates:
    adlib.adsorbate.xyz.make_xyz(ads, base_dir=xyz_dir)
