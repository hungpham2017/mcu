.. _plotting:

Plotting band structure
-----------------------
In python environment, band structure can be plotted by using **mcu.plot_band()** function

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_band()

To customize the band, one can modify these attributes:
* efermi: Fermi level. **Default**: from the vasprun.xml or OUTCAR (need to be specified in mcu.VASP()) 
* spin: Spin. For ISPIN = 1, spin = 0. For ISPIN = 2, spin = 1, 2. For LSORBIT = True, spin = 0, 1, 2, 3. **Default**: 0
* label: labels for high symmetric kpoints. For conventional band structure, label can be label = 'G-X-L-F'.
For hydrib functional, not only labels but also the coordinates need to be specified. **Default**: None. 
* save=False
* band_color=['#007acc','#808080','#808080']
* figsize=(6,6)
* figname='BAND'
* ylim=[-6,6]
* fontsize=18
* dpi=600
* format='png'
    

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