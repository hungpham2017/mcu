Welcome to mcu's documentation!
===============================
A package for post periodic wave function and crystallography analysis. mcu is designed for large scale analysis rather than for generating figures for one calculation.

A quick look
------------
A projected band structure can be plotted simply by:

.. code-block:: python
   :linenos:
   
   import mcu           
   mymcu = mcu.VASP()               # Define a mcu object
   mymcu.plot_pband(save=True)      # plot projected band structure, save=True to export an image

.. image:: ../image/MoS2.png
   :scale: 80 %
   :align: center

.. toctree::
    :glob:
    :maxdepth: 2
    :caption: Basic tutorials:
    :numbered:
   
    basic
    plotting
    *
   
Bugs and request new functions
------------------------------
* Bugs can be reported via openning an issue on my github or shooting me an email
* Suggest a new function? Shoot me an email

Contact
-------
Author: Hung Q. Pham
Email : pqh3.14@gmail.com