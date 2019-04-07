import mcu
run = mcu.VASP()
label = 'Y-G-R-X-G'
run.plot_pband(style=3,lm='pd',label=label,ylim=(-1.5,1.5),figsize=(8,6))