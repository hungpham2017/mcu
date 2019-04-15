import mcu
mymcu = mcu.VASP()

# Style = 1
mymcu.plot_dos(save=True, figname='Ni_horizontal', dpi=300)

# Style = 2 and spin = 'updown'
mymcu.plot_dos(spin = 'updown', style = 2, lm = ['Ni:s,dxy,dyz','Ni:p','Ni:dz2,dx2-y2'], save=True, figname='Ni_updown', dpi=300)