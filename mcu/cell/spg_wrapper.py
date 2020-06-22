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



def compare_cells(cell_1, cell_2, prec=1e-3):
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

def get_sym(cell, symprec=1e-3, verbose='short', angle_tolerance=-1.0, hall_number=0):
    '''Giving a spglib cell, return symmetry analysis
       This is wrapper for spglib.get_symmetry_dataset function 
    '''

    #Space group info
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=angle_tolerance, hall_number=hall_number)
    
    number = dataset['number'] 
    international = dataset['international']                    # equil. to international_short
    hall = dataset['hall']                                      # equil. to hall_symbol
    hall_number = dataset['hall_number']                        
    choice = dataset['choice'] 
    transformation_matrix = dataset['transformation_matrix'] 
    origin_shift = dataset['origin_shift'] 
    wyckoffs = dataset['wyckoffs'] 
    site_symmetry_symbols = dataset['site_symmetry_symbols']
    equivalent_atoms = dataset['equivalent_atoms']
    crystallographic_orbits = dataset['crystallographic_orbits']
    mapping_to_primitive = dataset['mapping_to_primitive']
    rotations = dataset['rotations']
    translations = dataset['translations']
    pointgroup = dataset['pointgroup']                          # equil. to pointgroup_international
    primitive_lattice = dataset['primitive_lattice']
    std_lattice = dataset['std_lattice']
    std_positions = dataset['std_positions']
    std_types = dataset['std_types']
    std_rotation_matrix = dataset['std_rotation_matrix']
    std_mapping_to_primitive =  dataset['std_mapping_to_primitive'] 
        
    # Get full summetry using the Hall number
    spacegroup_type = spglib.get_spacegroup_type(hall_number)
    number = spacegroup_type['number']
    international_short = spacegroup_type['international_short']
    international_full = spacegroup_type['international_full']
    international = spacegroup_type['international']
    schoenflies = spacegroup_type['schoenflies']
    hall_symbol = spacegroup_type['hall_symbol']
    pointgroup_schoenflies = spacegroup_type['pointgroup_schoenflies']
    pointgroup_international = spacegroup_type['pointgroup_international']
    arithmetic_crystal_class_number = spacegroup_type['arithmetic_crystal_class_number']
    arithmetic_crystal_class_symbol = spacegroup_type['arithmetic_crystal_class_symbol']

    atoms = utils.convert_atomtype(cell[2])
    coords = cell[1]
    std_cell = cell_to_std(dataset, message=False)
    primitive_cell = cell_to_primitive(dataset, message=False)
    
    if verbose == 'short':
        # Cell info
        if compare_cells(cell, std_cell) and compare_cells(cell, primitive_cell):
            misc.print_msg("This is a standard/primitive unit cell")
        elif compare_cells(cell, std_cell):
            misc.print_msg("This is a standard unit cell")
        elif compare_cells(cell, primitive_cell):
            misc.print_msg("This is a primitive cell") 
    
        # Basic info
        misc.print_msg("Space group  number          : {:d}".format(number))
        misc.print_msg("Short international symbol   : {:s}".format(international_short))
        misc.print_msg("Hall symbol                  : {:s}".format(hall_symbol))    
        misc.print_msg("Schoenflies symbol           : {:s}".format(schoenflies))
        misc.print_msg("International point group    : {:s}".format(pointgroup_international))
        misc.print_msg("Origin shift                 : {:4.3f} {:4.3f} {:4.3f}".format(origin_shift[0], origin_shift[1], origin_shift[2])) 

        # Atoms info
        # Atoms info
        misc.print_msg("===== Irreducible atoms list =====")
        misc.print_msg("  #  Atom        x          y          z")
        unique_atoms_idx = np.unique(equivalent_atoms, return_index=True)[1]
        for i, index in enumerate(unique_atoms_idx):
            coord = coords[index]
            misc.print_msg("{:>3d}  {:>3s}    {:>8.5f}   {:>8.5f}   {:>8.5f}".format(i+1, atoms[index], coord[0], coord[1], coord[2]))
        
    elif verbose == 'full':
        # Cell info
        if compare_cells(cell, std_cell) and compare_cells(cell, primitive_cell):
            misc.print_msg("This is a standard/primitive unit cell")
        elif compare_cells(cell, std_cell):
            misc.print_msg("This is a standard unit cell")
        elif compare_cells(cell, primitive_cell):
            misc.print_msg("This is a primitive cell") 
            
        # Basic info
        misc.print_msg("Space group  number          : {:d}".format(number))
        misc.print_msg("Short international symbol   : {:s}".format(international_short))
        misc.print_msg("Full international symbol    : {:s}".format(international_full))
        misc.print_msg("Hall number                  : {:d}".format(hall_number))
        misc.print_msg("Hall symbol                  : {:s}".format(hall_symbol))  
        misc.print_msg("Schoenflies symbol           : {:s}".format(schoenflies))
        misc.print_msg("Schoenflies point group      : {:s}".format(pointgroup_schoenflies))        
        misc.print_msg("International point group    : {:s}".format(pointgroup_international))
        misc.print_msg("Origin shift                 : {:4.3f} {:4.3f} {:4.3f}".format(origin_shift[0], origin_shift[1], origin_shift[2])) 
        
        # Atoms info
        misc.print_msg("===== Full atoms list =====")
        misc.print_msg("  #  Equil.  Atom        x         y         z     Wyckoffs   Site_symmetry")
        for i, atom in enumerate(atoms):
            misc.print_msg("{:>3d}  {:>3d}     {:>3s}    {:>8.5f}  {:>8.5f}  {:>8.5f}      {:>2s}          {:>5s}".format(i+1, equivalent_atoms[i] + 1, atoms[i], 
                            coords[i,0], coords[i,1], coords[i,2], wyckoffs[i], site_symmetry_symbols[i]))
    elif verbose is None:
        # Return an standard cell object with irreducible atoms
        dataset = spglib.get_symmetry_dataset(std_cell, symprec=symprec, angle_tolerance=angle_tolerance, hall_number=hall_number)
        equivalent_atoms = dataset['equivalent_atoms']
        irred_idx = np.unique(equivalent_atoms, return_index=True)[1]
        lattice = std_lattice
        irred_coord = std_positions[irred_idx,:]
        irred_label = std_types[irred_idx]
        irred_cell= (lattice, irred_coord, irred_label)  
        spacegroup = [number, international_short]
        
        # Get symmetry operators using the Hall number:
        # For some reason, the symmetry operators by dataset doesn't really match its space group number to write cif file
        symmetry = spglib.get_symmetry_from_database(hall_number)
        rotations = symmetry['rotations'] 
        translations = symmetry['translations'] 
    
        return spacegroup, irred_cell, rotations, translations
    
def cell_to_std(cell_or_dataset, symprec=1e-3, angle_tolerance=-1.0, hall_number=0, message=True):
    '''Convert a cell to a standard cell'''

    if isinstance(cell_or_dataset, dict):
        dataset = cell_or_dataset
    else:
        cell = tuple(cell_or_dataset)
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=angle_tolerance, hall_number=hall_number)
        
    # Construct a std cell object    
    std_lattice = dataset['std_lattice']
    std_positions = dataset['std_positions']
    std_types = dataset['std_types']
    std_cell = (std_lattice, std_positions, std_types)
    
    if message:
        if compare_cells(cell, std_cell):
            print('The unit cell was a standard cell. However, the standard cell computed by spglib is returned, it is maybe not the same as the provided unit cell')
        else:
            print('The unit cell was transformed to a standard cell')  
        
    return std_cell

def cell_to_primitive(cell_or_dataset, symprec=1e-3, angle_tolerance=-1.0, hall_number=0, message=True):
    '''Convert a cell to a primitive cell'''
    
    if isinstance(cell_or_dataset, dict):
        dataset = cell_or_dataset
    else:
        cell = tuple(cell_or_dataset)
        dataset = spglib.get_symmetry_dataset(cell, symprec=symprec, angle_tolerance=angle_tolerance, hall_number=hall_number)
        
    # Construct a primitive cell object    
    std_lattice = dataset['std_lattice']
    std_positions = dataset['std_positions']
    std_types = dataset['std_types']
    std_cell = (std_lattice, std_positions, std_types)

    lattice, scaled_positions, numbers = spglib.standardize_cell(std_cell, to_primitive=True, no_idealize=True, symprec=symprec)
    primitive_cell = (lattice, scaled_positions, numbers)
    
    if message:
        if compare_cells(cell, primitive_cell):
            print('The unit cell was a primitive cell. However, the primitive cell computed by spglib is returned, it is maybe not the same as the provided unit cell')
        else:
            print('The unit cell was transformed to a primitive cell') 
            
    return primitive_cell   
    
get_symmetry_from_database = spglib.get_symmetry_from_database
get_ir_reciprocal_mesh = spglib.get_ir_reciprocal_mesh