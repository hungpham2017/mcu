__version__ = '1.0.0'
from . import vasp


VASP = vasp.main.VASP
LOCPOT = vasp.main.LOCPOT
POSCAR = vasp.main.POSCAR
vasprun = vasp.read.vasprun
OUTCAR = vasp.read.OUTCAR
make_KPOINTS = vasp.utils.get_1Dkpath
read_WAVEDER = vasp.utils.read_WAVEDER
read_WAVEDERF = vasp.utils.read_WAVEDERF

