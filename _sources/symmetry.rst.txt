.. _symmetry:

Spacegroup analysis 
===================

**mcu** can analyze spacegroup symmetry of the crystal. 
You can analysze whether the crystal is a standard unit cell or a primitive unit cell, then transform it between them.
**mcu** can export variety of crystal structure formats, such as cif, POSCAR, xsf. For cif format, you can either
get the structure in its spacegroup or in P1 symmetry. These analyses can be started with a **mcu.VASP** or **mcu.POSCAR** object.
Here is an example:

.. code-block:: python
    :linenos:
   
    import mcu
    
    # Define a POSCAR object with POSCAR/CONTCAR file    
    mymcu = mcu.POSCAR()
    
    # Analyze spacegroup
    mymcu.get_symmetry()
    
    # Transformto a primitive unit cell
    mymcu.to_primcell()
    mymcu.get_symmetry()  

    # Transformto a standard unit cell
    mymcu.to_stdcell()
    mymcu.get_symmetry()      
    mymcu.write_poscar()
    
    # Export to cif. The current symmetry if analyzed will be included in cif file.
    # if get_symmetry() is not executed yet, then a P1 structure will be returned.
    mymcu.write_cif()
    
    # Export a P1 structure even though the spacegroup was analyzed
    mymcu.write_cif(symmetry=False)    

By default, the **mymcu.cell** (read from POSCAR/CONTCAR or vasprun.xml) is used for the above analyses.
The **mymcu.cell** is in **spglib** format so you actually can pass it to the **spglib** module to do other analyses.
Moreover, you can analyze symmtry for any structure of interest.

.. code-block:: python
    :linenos:
   
    import mcu
    
    # Define a POSCAR object with POSCAR/CONTCAR file    
    mymcu = mcu.POSCAR()
    
    cell = mymcu.cell       # in general, any cell type
    mymcu.get_symmetry(cell)

    # mymcu.cell is transformaed to a primitive cell    
    mymcu.to_primcell()
    prim_cell = mymcu.cell
    mymcu.get_symmetry(prim_cell)

    # export a cif file for a cell
    mymcu.write_cif(prim_cell)      

You can read any cif (e.g. cif from XRD refinement) from mcu and export it to your desire format (POSCAR or xsf).
One should be careful about the possible disosder (e.g. very small distance between atoms) from the XRD cif file.
 
.. code-block:: python
    :linenos:
   
    import mcu
    mycif = mcu.CIF('yourcif.cif')      #if a file is not provided, then the first cif file (if you run command ls) will be processed 
    mycif.write_poscar()
    
