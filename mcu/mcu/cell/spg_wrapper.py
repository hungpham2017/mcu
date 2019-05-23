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
import spglib
from mcu.utils import misc
from mcu.cell import parameters, utils



def compare_cells(cell1, cell2, prec=1e-5):
    '''Compare two spglib cell if they are the same'''
    lattice_diff = np.linalg.norm(np.asarray(cell1[0]) - np.asarray(cell2[0]))
    
    same = False
    if lattice_diff < prec:
        pos1 = np.asarray(cell1[1])
        pos2 = np.asarray(cell2[1]) 
        if pos1.shape[0] == pos2.shape[0]:
            if len(cell1[2]) == len(cell2[2]):
                atom_type = (cell1[2] == cell2[2]).all()
                if atom_type == True:
                    sym1 = spglib.get_symmetry(cell1, prec)
                    sym2 = spglib.get_symmetry(cell2, prec) 
                    rota1 = np.asarray(sym1['rotations'])
                    trans1 = np.asarray(sym1['translations'])
                    rota2 = np.asarray(sym2['rotations']) 
                    trans2 = np.asarray(sym2['translations'])
                    equi1 = np.asarray(sym1['equivalent_atoms'])
                    equi2 = np.asarray(sym2['equivalent_atoms'])                    
                    if rota1.shape[0] == rota2.shape[0]: 
                        diff_rotation = np.linalg.norm(rota1 - rota2) < prec
                        if diff_rotation and (trans1.shape[0] == trans2.shape[0]):
                            diff_translation = np.linalg.norm(trans1 - trans2) < prec
                            if diff_translation and (equi1.shape[0] == equi2.shape[0]):
                                diff_equi = np.linalg.norm(equi1 - equi2) < prec
                                if diff_equi:
                                    same = True
                                
    return same

def get_sym(cell, symprec=1e-5, print_atom=False, export_operator=False):
    
    #Space group info
    intl_label, number = spglib.get_spacegroup(cell, symprec).split(' ')
    number = number.split("(")[1].split(")")[0]
    Schoenflies_label = spglib.get_spacegroup(cell, symprec, symbol_type=1).split(' ')[0]
    sym = spglib.get_symmetry(cell, symprec)
    rotations = sym['rotations'] 
    translations = sym['translations'] 
    equi_atoms = sym['equivalent_atoms'] 
    
    is_std = is_prim = False
    std_cell = spglib.refine_cell(cell, symprec)
    prim_cell = spglib.find_primitive(cell, symprec)    
    if compare_cells(cell, std_cell): is_std = True
    if compare_cells(cell, prim_cell): is_prim = True
    
    atoms = utils.convert_atomtype(cell[2])
    if export_operator == False:
        #Cell info
        std_cell = spglib.refine_cell(cell, symprec)
        prim_cell = spglib.find_primitive(cell, symprec)
        #Print
        misc.print_msg("Spacegroup  number           : %s" % (number))
        misc.print_msg("Short Internatmiscnal symbol : %s" % (intl_label))
        misc.print_msg("Schoenflies symbol           : %s" % (Schoenflies_label))
        if print_atom == True:
            misc.print_msg("Atoms list (No.   Sym   symbol):")
            for i, atom in enumerate(atoms):
                misc.print_msg("%3d   %3d   %s" %(i+1, equi_atoms[i]+1, atom))
            misc.print_msg("Irreducible atoms:")
            for i, index in enumerate(np.unique(equi_atoms)):
                misc.print_msg("%3d  %3s    %7f5 %7f5 %7f5" %(i+1, atoms[index], cell[1][index][0], cell[1][index][1], cell[1][index][2]))
        misc.print_msg("Number of irreducible atoms  : %d" % (np.unique(equi_atoms).shape[0]))
        misc.print_msg("Standard cell                : %r" % (is_std))
        misc.print_msg("Primitive cell               : %r" % (is_prim)) 
        return is_std, is_prim
    else:
        return (number, intl_label), equi_atoms, rotations, translations

def cell_to_std(cell, symprec=1e-5):
    std_cell = spglib.refine_cell(cell, symprec) 
    if compare_cells(cell, std_cell):
        print('Unit cell is already a standard cell. Nothing changes')
        return cell
    else:
        print('Unit cell was transformed to a standard cell')  
        return std_cell

def cell_to_prim(cell, symprec=1e-5):
    prim_cell = spglib.find_primitive(cell, symprec)   
    if compare_cells(cell, prim_cell):
        print('Unit cell is already a primitive cell. Nothing changes') 
        return cell        
    else:
        print('Unit cell was transformed to a primitive cell')
        return prim_cell

get_symmetry_from_database = spglib.get_symmetry_from_database       

    
    
                      
        

        