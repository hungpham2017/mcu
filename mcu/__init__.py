#!/usr/bin/env python

__version__ = '0.2'
from . import vasp, cell, wannier90, cp2k, crystal

# Quick shortcuts to VASP tools
VASP = vasp.vasprun.main
LOCPOT = vasp.locpot.main
POSCAR = vasp.poscar.main
WAVECAR = vasp.wavecar.main
CIF = cell.cell_io.cif
vasprun = vasp.vasp_io.vasprun
OUTCAR = vasp.vasp_io.OUTCAR
make_KPOINTS = vasp.utils.get_1Dkpath
read_WAVEDER = vasp.utils.read_WAVEDER
read_WAVEDERF = vasp.utils.read_WAVEDERF

CELL = cell.main.CELL

# Quick shortcuts to wannier90 tools
W90 = wannier90.w90.main
read_unk = wannier90.utils.read_unk
read_U_matrix = wannier90.utils.read_U_matrix

# Quick shortcuts to cp2k tools
CP2K = cp2k.cp2k.main

# Quick shortcuts to cp2k tools
CRYSTAL = crystal.crystal.main
