#!/usr/bin/env python

__version__ = '0.9.0'
from . import utils
from . import cell
from . import vasp, cp2k, crystal, elk, pscf, qe, wannier90

# Quick shortcuts to VASP tools
VASP = vasp.vasprun.main

# TODO: The following syntax looks confusing and messing, make it simpler!!!
LOCPOT = vasp.locpot.main
POSCAR = vasp.poscar.main
WAVECAR = vasp.wavecar.main
vasprun = vasp.vasp_io.XML
make_KPOINTS = vasp.utils.get_1Dkpath
read_WAVEDER = vasp.utils.read_WAVEDER
read_WAVEDERF = vasp.utils.read_WAVEDERF


# Quick shortcuts to cell tools
CELL = cell.cell.main

# Quick shortcuts to wannier90 tools
W90 = wannier90.w90.main
read_unk = wannier90.utils.read_unk
read_U_matrix = wannier90.utils.read_U_matrix

# Quick shortcuts to cp2k tools
CP2K = cp2k.cp2k.main

# Quick shortcuts to cp2k tools
CRYSTAL = crystal.crystal.main

# Quick shortcuts to QE tools
QE = qe.qe.main

# Quick shortcuts to QE tools
PYSCF = pscf.pscf.main

# Quick shortcuts to QE tools
ELK = elk.elk.main
