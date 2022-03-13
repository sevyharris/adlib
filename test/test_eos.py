import adlib.bulk.eos


bulk_dir = '/work/westgroup/harris.se/espresso/adlib/test/bulk'

adlib.bulk.eos.setup_eos(bulk_dir, metal='Cu', lattice_constant_guess=3.6)
adlib.bulk.eos.run_eos(bulk_dir)

