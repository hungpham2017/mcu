.. _basic:

Basic functions in mcu 
===============================

**Compute bandgap and Fermi level:**

.. code-block:: python
   :linenos:
   
   import mcu
   mymcu = mcu.VASP()    # Define an mcu object
   mymcu.get_bandgap()   # Get the bandgap
   mymcu.efermi          # Get Fermi level                         

**WAVEDER and WAVEDERF file can be read from mcu:**

.. code-block:: python
   :linenos:
   
   import mcu
   cder, nodesn_i_dielectric_function, wplasmon = mcu.read_WAVEDER(waveder = 'WAVEDER')
   cder = mcu.read_WAVEDERF(wavederf = 'WAVEDERF')   