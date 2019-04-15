.. _plotting:

Plotting band structure
-----------------------
In python environment, band structure can be plotted by calling the **mcu.plot_band()** function

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_band()

To customize the band, one can modify some of these attributes:

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               
   mymcu.plot_band(spin=0, save=True, figsize=(6,6), dpi=600, format='png')
    
All parameters and their defaults of **plot_band** function are given below.

Parameters
++++++++++
efermi : float
         Fermi level.

    
    
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