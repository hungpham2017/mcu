.. _basic:

Basic functions in mcu 
===============================
#. Compute bandgap:

.. code-block:: python
   :linenos:
   
   import mcu
   mymcu = mcu.VASP()    # Define an mcu object
   mymcu.get_bandgap()   # Get the bandgap
   mymcu.efermi          # Get Fermi level                         


   