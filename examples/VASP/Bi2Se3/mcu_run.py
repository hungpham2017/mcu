import mcu

mymcu = mcu.VASP()
klabel = 'G-Z-F-G-L'

# Plot band structure:
mymcu.plot_pband(klabel = klabel)

# Plot projected band structure:
mymcu.plot_pband(lm='p;s', legend=['p','s'], klabel=klabel, color=['#00ccff','#ff0000'],alpha=0.4,ylim=(-2,2),legend_size=1.2)

# Style 3
mymcu.plot_pband(lm='s',color='#ff0000',alpha=0.4,label=label,ylim=(-2,2),figsize=(8,6),
scale=10,legend='s',legend_size=1.2,marker='h', save=True, figname='Bi2Se3_b', dpi=300)