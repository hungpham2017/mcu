import mcu
run = mcu.VASP()
label = 'Y-G-R-X-G'

# Style 2
run.plot_pband(style=2,lm=['Mo:d','S:p'],color=['#00ccff','#ff0000'],alpha=0.4,label=label,ylim=(-1.5,1.5),figsize=(8,6))
# Style 3
run.plot_pband(style=3,lm='pd',label=label,ylim=(-1.5,1.5),figsize=(8,6))