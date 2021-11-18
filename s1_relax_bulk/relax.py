from ase.build import bulk, fcc111, add_adsorbate
from ase.io import write, read
from ase import Atoms
from ase.optimize import LBFGS
from ase.calculators.espresso import Espresso
from ase.io.espresso import write_espresso_in


out_file = 'bulk.xyz'
cu_bulk = bulk('Cu', crystalstructure='fcc', a=3.0, cubic=True)
print(cu_bulk.get_volume())
print(cu_bulk.cell.get_bravais_lattice())

write(out_file, cu_bulk)
pseudopotentials = {'Cu': 'Cu.pbe-dn-kjpaw_psl.1.0.0.UPF'}

input_settings = {
    'control': {
        'calculation': 'vc-relax',
        'forc_conv_thr': 0.001
    },
    'system': {
        'occupations': 'smearing',
        'degauss': 0.1,
        # 'A': 2.468,
        # 'C': 8.685,
        # 'nat': 4,
        # 'ntyp': 1,
        'ecutwfc': 25,
        'ecutrho': 400
    },
    'ions': {
        'ion_dynamics': 'bfgs'
    },
    'cell': {
        'cell_dynamics': 'bfgs',
        'press': 0.0,
        'press_conv_thr': 0.5,
    }

}


calc = Espresso(pseudopotentials=pseudopotentials,
                tstress=True, tprnfor=True, kpts=(1, 1, 1),
                pseudo_dir='.',
                input_data=input_settings)
cu_bulk.calc = calc

# ucf = UnitCellFilter(rocksalt)
opt = LBFGS(cu_bulk)
opt.run(fmax=0.005)

# compute the optimized cubic lattice constant
print(cu_bulk.cell.get_bravais_lattice())
print(cu_bulk.get_volume())

write('cu_opt.xyz', cu_bulk)

