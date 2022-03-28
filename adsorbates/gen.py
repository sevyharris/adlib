from ase.io import write
from ase import Atoms
from ase.collections import g2


def gen_CO2():
    # CO2 - need to make sure the orientation is relative to the z normal of the slab
    # d = 1.0
    # co2 = Atoms('CO2', positions=[(0, 0, 0), (d, 0, 0), ((-d, 0, 0))])
    co2 = g2['CO2']
    co2.rotate(90, 'y')
    write('CO2.xyz', co2)


def gen_CO():
    # d = 1.0
    # co = Atoms('CO', positions=[(0, 0, 0), (d, 0, 0)])
    co = g2['CO']
    co.rotate(90, 'y')
    write('CO.xyz', co)


def gen_O2():
    # d = 1.0
    # o2 = Atoms('O2', positions=[(0, 0, 0), (d, 0, 0)])
    o2 = g2['O2']
    o2.rotate(90, 'y')
    write('O2.xyz', o2)


def gen_O():
    o = Atoms('O')
    write('O.xyz', o)


def gen_H2():
    # d = 1.0
    # h2 = Atoms('H2', positions=[(0, 0, 0), (d, 0, 0)])
    h2 = g2['H2']
    h2.rotate(90, 'y')
    write('H2.xyz', h2)


def gen_H2O():
    d = 1.0
    # h2o = Atoms('H2O', positions=[(d, 0, 0), (0, 0, 0), ((-d, 0, 0))])
    h2o = g2['H2O']
    h2o.rotate(90, 'y')
    write('H2O.xyz', h2o)


def gen_H():
    h = Atoms('H')
    write('H.xyz', h)


def gen_OH():
    # d = 1.0
    # oh = Atoms('OH', positions=[(0, 0, 0), (d, 0, 0)])
    oh = g2['OH']
    oh.rotate(90, 'y')
    write('OH.xyz', oh)


def gen_N2():
    # d = 1.0
    # n2 = Atoms('N2', positions=[(0, 0, 0), (d, 0, 0)])
    n2 = g2['N2']
    n2.rotate(90, 'y')
    write('N2.xyz', n2)


def gen_CH4():
    # ch4 = Atoms('CH4')
    ch4 = g2['CH4']
    write('CH4.xyz', ch4)


def gen_CH():
    ch = g2['CH']
    write('CH.xyz', ch)


def gen_CH3():
    ch3 = g2['CH3']
    write('CH3.xyz', ch3)



if __name__ == '__main__':
    #gen_CO2()
    #gen_CO()
    #gen_O2()
    #gen_O()
    #gen_H2()
    #gen_H2O()
    #gen_H()
    #gen_OH()
    #gen_N2()
    #gen_CH4()
    gen_CH()
    gen_CH3()
