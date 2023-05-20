import mcu

mymcu = mcu.VASP()

def test_plotband():
    '''Testing plot_band function'''
    mymcu.get_bandgap()
    out = mymcu.plot_pband(save=True)
    
def test_sym():
    mymcu.get_symmetry()
