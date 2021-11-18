Requires:
Quantum Espresso
ASE

Steps:
1. vc-relax a bulk unit cell to get the bulk lattice constant for CU using that functional (PBE to start)
2. Construct a clean slab using the results from 1. Fix the bottom layers (everything but the top two layers) and relax the slab
3. Relax isolated adsorbate (CO2) - do I need same dimensions as unit cell for slab? Probably not, but could make combining things easier
4. Find equilibrium distance of CO2 from top site keeping rigid geometries of slab and bulk
5. Relax whole system


