__version__ = '1.0.0'
from . import vasp


VASP = vasp.main.VASP
vasprun = vasp.vasprun
read_WAVEDER = vasp.utils.read_WAVEDER
read_WAVEDERF = vasp.utils.read_WAVEDERF
