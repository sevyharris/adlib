import adlib.bulk.vcrelax


bulk_dir = '/work/westgroup/harris.se/espresso/adlib/test/bulk'

adlib.bulk.vcrelax.setup_vc_relax(bulk_dir, metal='Cu', lattice_constant_guess=3.6)
adlib.bulk.vcrelax.run_vc_relax(bulk_dir)
