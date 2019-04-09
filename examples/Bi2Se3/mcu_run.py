import mcu
run = mcu.VASP()
label = 'G-Z-F-G-L'

# Style 2
run.plot_pband(style=1,lm='ps',color=['#00ccff','#ff0000'],alpha=0.4,label=label,ylim=(-2,2),figsize=(8,6),legend=['p','s'],legend_size=1.2)
# Style 3
run.plot_pband(style=1,lm='s',color='#ff0000',alpha=0.4,label=label,ylim=(-2,2),figsize=(8,6),scale=10,legend='s',legend_size=1.2,marker='h')