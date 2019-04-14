.. _plotting:

Plotting band structure
-----------------------
In python environment, band structure can be plotted by calling the **mcu.plot_band()** function

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_band()

To customize the band, one can modify these attributes:

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_band(efermi=None, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'], figsize=(6,6), figname='BAND', ylim=[-6,6], fontsize=18, dpi=600, format='png')
    

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