import adlib.bulk.convergence


bulk_dir = '/work/westgroup/harris.se/espresso/adlib/test/bulk'
lattice_constant = 3.57
adlib.bulk.convergence.setup_ecutwfc_converge(bulk_dir, lattice_constant)
adlib.bulk.convergence.setup_kpts_converge(bulk_dir, lattice_constant)
adlib.bulk.convergence.setup_smear_converge(bulk_dir, lattice_constant)
