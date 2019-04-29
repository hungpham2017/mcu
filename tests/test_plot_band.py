import mcu

mymcu = mcu.VASP()

def test_plotband():
    out = mymcu.plot_pband(save=True)
