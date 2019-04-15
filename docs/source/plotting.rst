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
    
All the attributes and their defaults are given below.

def plot_band(efermi=None, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'],
                    figsize=(6,6), figname='BAND', ylim=[-6,6], fontsize=18, dpi=600, format='png'):
    """Plot band structure
    
    :param efermi: None
    :param spin: 0
    :param label: None
    :param save: False
    :param band_color: ['#007acc','#808080','#808080']
    :param figsize: (6,6)
    :param figname: 'BAND'
    :param ylim: [-6,6]
    :param fontsize: 18
    :param dpi: 600
    :param format: 'png'
    """
    
    
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