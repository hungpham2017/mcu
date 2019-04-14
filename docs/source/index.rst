Welcome to mcu's documentation!
===============================
A package for post periodic wave function and crystallography analysis. mcu is designed for large scale analysis rather than for generating figures for one calculation.

A quick look
============
A projected band structure can be plotted simply by:

.. code-block:: python
   :linenos:
   
   import mcu
   mymcu = mcu.VASP()
   mymcu.plot_pband()

.. image:: ../image/MoS2.png

.. toctree::
   :maxdepth: 2
   :caption: Basic tutorials:
   
   basic
   plotting
   