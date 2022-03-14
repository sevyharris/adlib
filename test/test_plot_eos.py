import adlib.bulk.eos


bulk_dir = '/work/westgroup/harris.se/espresso/adlib/test/bulk'

lattice_constant = adlib.bulk.eos.analyze_eos(bulk_dir)
print(f'Lattice constant: {lattice_constant}')

adlib.bulk.eos.plot_eos(bulk_dir)
