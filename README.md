# Description
adlib (Adsorption Library) is a collection of scripts used to generate [Quantum Espresso](https://www.quantum-espresso.org/) jobs for computing adsorption energies

# How to Install
Download the source code:

`git clone https://github.com/sevyharris/adlib.git`

Install with pip:

`pip install -e adlib`

# Brief Overview
adlib contains functions to generate python scripts that use [ase](https://wiki.fysik.dtu.dk/ase/) and [Quantum Espresso](https://www.quantum-espresso.org/) to compute energies for adsorption.

These functions are divided into
* bulk (3D repeating metal)
* slab (2D repeating metal)
* adsorbate (isolated gas molecule)
* system (combined slab and adsorbate).


adlib follows this workflow to compute binding energies:
1. Compute bulk lattice constant (do vc relax for initial guess, then eos coarse and eos fine to narrow in on minimum energy)
2. Optimize slab geometry and compute energy
3. Optimize adsorbate geometry and compute energy
4. Place adsorbate on slab, optimize geometry, compute energy

The binding energy is then:

system - (slab + adsorbate)


adlib assumes the results will be organized in a particular directory structure, detailed in the comments for each function
<!--- TODO - describe this in more detail here --->

If you just want some example files to look at, see the [basepys](https://github.com/sevyharris/adlib/tree/main/basepys) folder

# Documentation
See the documentation for more instructions and examples
[https://sevyharris.github.io/adlib/](https://sevyharris.github.io/adlib/)
