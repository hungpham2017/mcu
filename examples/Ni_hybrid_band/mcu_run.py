import mcu
run = mcu.VASP()

# Get E_fermi
run.vasprun.get_dos()
efermi = run.vasprun.efermi

#Get band gap
run.get_bandgap(efermi)

#Plot band structure
label = [['G',0.0,0.0,0.00],['X',0.250000,0.500000,0.250000]]
run.plot_band(efermi,label=label,hybridXC=True)