Welcome to mcu's documentation!
===============================
A package for periodic wavefunction and crystallography analysis. 
**mcu** is designed to support large scale analysis and topological descriptions for periodic wavefunction.

A quick look
------------
A projected band structure can be plotted simply by:

.. code-block:: python
    :linenos:
   
    import mcu
    mymcu = mcu.VASP()
    mymcu.plot_pband()

.. image:: ../image/MoS2.png
    :scale: 60 %
    :align: center
    
**mcu** faciliates the setup and plotting for 2D band structure. 
For example, one can visualize the Dirac cones of graphene:

.. code-block:: python
    :linenos:
   
    import mcu
    mymcu = mcu.VASP()
    mymcu.plot_band2D()

.. image:: ../image/gra_band2D.png
    :scale: 60 %
    :align: center
    
or the spin texture:

.. image:: ../image/elecpot.png
    :scale: 20 %
    :align: center

    
Content:
--------

.. toctree::
    :maxdepth: 2

    feature
    install
    tutorial
    gallery
    contact

    
    
