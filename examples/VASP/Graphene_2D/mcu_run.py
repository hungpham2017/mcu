import mcu

# Generate 2D k-mesh from POSCAR
run1 = mcu.POSCAR()
run1.get_2D_kmesh(origin=[0,0,0], krange=[0.2,0.2], plane='xy', npoint=[21,21])

# Plot band structure
run2 = mcu.VASP()
run2.plot_band2D()
