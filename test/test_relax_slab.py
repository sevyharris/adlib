import adlib.slab.calc


slab_dir = '/work/westgroup/harris.se/espresso/adlib/test/slab'

adlib.slab.calc.make_relax_script(slab_dir, 3.6)
adlib.slab.calc.make_run_relax_script(slab_dir)
adlib.slab.calc.run_relax_slab(slab_dir)
