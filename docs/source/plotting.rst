.. _plotting:

Plotting band structure
-----------------------
In python environment, band structure can be plotted by calling the **mcu.plot_band()** function

.. code-block:: python
    :linenos:
   
    import mcu           
    mymcu = mcu.VASP()    # or mymcu = mcu.VASP(path='path-to-vasprun', vaspruns='vasprun')             
    mymcu.plot_band()

To customize the band, one can modify some of these attributes. For mcu/exampple/MoS2, you can run:

.. code-block:: python
    :linenos:
   
import mcu           
mymcu = mcu.VASP()   
mymcu.plot_band(spin=0, save=True, label='Y-G-R-X-G', ylim=(-3,3), figsize=(3,3), dpi=300, format='png')
    
You should get:

.. image:: ../image/plot_band_MoS2_a.png
    :scale: 80 %
    :align: center
    
All parameters and their defaults of **plot_band** function are given below. Most of the parameters are passed to matplotlib functions.
So more information can be found in matplotlib docs.

Parameters
~~~~~~~~~~
efermi : float
    * Default: fermi level from vasprun.xml or OUTCAR
    * User can shift the Fermi level to a value
spin : int
    * Default: 0
    * If ISPIN = 1: spin = 0
    * If ISPIN = 2: spin = 0 (Up spin) or 1 (Down spin)
label : str or list
    * Default: None 
    * For conventional band structure, e.g. label = 'X-G-Y-L-G'
    * For hydrib functional band structure, e.g. label = [['L',0.50,0.50,0.50],['G',0.0,0.0,0.00],['X',0.5,0.0,0.50],['W',0.50,0.25,0.75]]
save : bool 
    * Default: False
    * True to save to an image    
band_color: list
    * Default: ['#007acc','#808080','#808080']
    * Three color codes indice color of band curve, kpoint grid, Fermi level, respectively. 
    * Exp: ['k','#808080','r']
figsize : tuple or list
    * Default: De(6,6)
    * Size of image in inch
figname : str
    * Default: 'BAND'
    * Name of the image
ylim : list
    * Default: [-6,6]
    * Limit range for energy axis in eV
fontsize : int
    * Default: 18
    * Font size
dpi : int
    * Default: 600
    * Resolution of the image 
format : str
    * Default: 'png'
    * Extension of the image
    
    
Plotting projected band structure
---------------------------------
In python environment, band structure can be plotted by calling the **mcu.plot_band()** function

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_pband()
   

Plotting density of states
--------------------------
For DOS, the total DOS is always shown together with projected DOS (if computed).

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_dos()

.. toctree::
    :glob:
    :maxdepth: 2
    :caption: Contents: