#!/usr/bin/env python
'''
mcu: Modeling and Crystallographic Utilities
Copyright (C) 2019 Hung Q. Pham. All Rights Reserved.

Licensed under the Apache License, Versmiscn 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITmiscNS OF ANY KIND, either express or implied.
See the License for the specific language governing permissmiscns and
limitatmiscns under the License.

Email: Hung Q. Pham <pqh3.14@gmail.com>
'''

import numpy as np
import spglib, copy
from ..utils import misc
from . import parameters, utils



def compare_cells(cell_1, cell_2, prec=1e-5):
    '''Compare two spglib cell if they are the same'''
    
    # spglib only work with tuple
    lattice_1, coord_1, atom_1 = tuple(cell_1)
    lattice_2, coord_2, atom_2 = tuple(cell_2)     
    is_the_same = True
    
    # Check volume of the unit cell
    volume_1 = np.linalg.det(lattice_1)
    volume_2 = np.linalg.det(lattice_2)
    if abs(volume_1 - volume_2) > prec:
        is_the_same = False
    if len(atom_1) != len(atom_2):
        is_the_same = False
            
    return is_the_same

def get_sym(cell, symprec=1e-5, print_atom=False, print_analysis=True):
    '''Giving a spglib cell, return symmetry analysis'''
    
    # spglib only work with tuple
    cell = tuple(cell)
    
    #Space group info
    spg_label, spg_number = spglib.get_spacegroup(cell, symprec).split(' ')
    spg_number = spg_number.split("(")[1].split(")")[0]
    Schoenflies_label = spglib.get_spacegroup(cell, symprec, symbol_type=1).split(' ')[0]
    sym = spglib.get_symmetry(cell, symprec)
    rotations = sym['rotations'] 
    translations = sym['translations'] 
    equi_atoms = sym['equivalent_atoms']  
    atoms = utils.convert_atomtype(cell[2])
    
    if print_analysis is True:
        #Cell info
        std_cell = spglib.refine_cell(cell, symprec)
        prim_cell = spglib.find_primitive(cell, symprec)
        if compare_cells(cell, std_cell) and compare_cells(cell, prim_cell):
            misc.print_msg("This is a standard/primitive unit cell")
        elif compare_cells(cell, std_cell):
            misc.print_msg("This is a standard unit cell")
        elif compare_cells(cell, prim_cell):
            misc.print_msg("This is a primitive cell") 
            
            
        #Print
        misc.print_msg("Spacegroup  number           : %s" % (spg_number))
        misc.print_msg("Short International symbol   : %s" % (spg_label))
        misc.print_msg("Schoenflies symbol           : %s" % (Schoenflies_label))
        if print_atom is True:
            misc.print_msg("Atoms list (No. - Sym - Symbol):")
            for i, atom in enumerate(atoms):
                misc.print_msg("%3d   %3d   %s" %(i+1, equi_atoms[i]+1, atom))
            misc.print_msg("Irreducible atoms:")
            for i, index in enumerate(np.unique(equi_atoms)):
                coord = cell[1][index]
                misc.print_msg("%3d  %3s    %7f5 %7f5 %7f5" %(i+1, atoms[index], coord[0], coord[1], coord[2]))
        misc.print_msg("Number of irreducible atoms  : %d" % (np.unique(equi_atoms).shape[0]))
    else:
        # Return an standard cell object with irreducible atoms
        
        irred_idx = np.unique(equi_atoms)
        lattice = cell[0]
        irred_coord = cell[1][irred_idx,:]
        irred_label = np.array(cell[2])[irred_idx]
        irred_cell= (lattice, irred_coord, irred_label)  
        spacegroup = [int(spg_number), spg_label]
        return spacegroup, irred_cell, rotations, translations
    
def cell_to_std(cell, symprec=1e-5):
    '''Convert a cell to a standard cell'''
    
    # spglib only work with tuple
    cell = tuple(cell)
    std_cell = spglib.refine_cell(cell, symprec) 
    if compare_cells(cell, std_cell):
        print('The unit cell was a standard cell. However, the standard cell computed by spglib is returned, it is maybe not the same as the provided unit cell')
    else:
        print('The unit cell was transformed to a standard cell')  
    return std_cell

def cell_to_prim(cell, symprec=1e-5):
    '''Convert a cell to a primitive cell'''
    
    # spglib only work with tuple
    cell = tuple(cell)
    prim_cell = spglib.find_primitive(cell, symprec) 
    if compare_cells(cell, prim_cell):
        print('The unit cell was a primitive cell. However, the primitive cell computed by spglib is returned, it is maybe not the same as the provided unit cell')
    else:
        print('The unit cell was transformed to a primitive cell') 
    return prim_cell

get_symmetry_from_database = spglib.get_symmetry_from_database       

    
    
                      
        

        