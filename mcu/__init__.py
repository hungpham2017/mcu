#!/usr/bin/env python

__version__ = '0.0.1'
from . import vasp
from . import cell

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