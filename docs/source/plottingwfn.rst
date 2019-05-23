.. _plottingwfn:

..
    ///////////////Plotting the U_n^k(r)///////////////   
    
Plotting the periodic part unk of a Bloch state
===============================================

POSCAR is needed to to get the lattice information. 
the following command generates a k-mesh (-0.2 <= x <= 0.2 and -0.2 <= y <= 0.2) on the xy plan and around the :math:`\Gamma` point with 21 points along each dimension.

.. code-block:: python 
    :linenos:
   
    import mcu           
    mymcu = mcu.WAVECAR() 
    unk = mymcu.get_unk()       # Generate the unk
    mymcu.write_vesta(unk)      # Save to a VESTA file for visualization



Parameters for **get_unk**
~~~~~~~~~~~~~~~~~~~~~~~~~~    
spin : int
    * Default: 0
    * If ISPIN = 1 or LSORBIT = True: spin = 0
    * If ISPIN = 2: spin = 0 (up spin) or 1 (down spin)
kpt : int
    * Default: 1
    * The first k-point in the k-point list
band : 1
    * Default: 1
    * The first band  
ngrid : list
    * Default: 2*ngrid with ngrid is the estimated minimum planewave
    * The mesh used in FFT of planewave coefficients to get the unk
norm_u : bool
    * Default: True
    * Normalize the unk by 1/sqrt(N_G) if N_G if the number of grid point
norm_c: bool
    * Default: False
    * Whether to normalize the planewave coefficients (C_G) such that <C_G|C_G> = 1 
    
Parameters for **write_vesta**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
unk : unk object
    * Default: have to provide one unk to plot
realonly: bool
    * Default: False
    * Generate only real part of the unk
    
poscar: POSCAR file
    * Default: POSCAR in the current directory
filename : str
    * Default: 'unk'
    * The unk_r.vasp (real) and unk_i.vasp (imaginary) will be generatated
ncol : int 
    * Default: 10
    * The number of column of the unk values in the VESTA file