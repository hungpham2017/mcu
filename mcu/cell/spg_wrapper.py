#!/usr/bin/env python
'''
mcu: Modeling and Crystallographic Utilities
Copyright (C) 2019 Hung Q. Pham. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Email: Hung Q. Pham <pqh3.14@gmail.com>
'''

import numpy as np
import spglib

        
class SPG:
    def __init__(self, cell, symprec=1e-5):
        self.cell = cell
        self.symprec = symprec
    
    def get_error(self):    
        return spglib.get_error_message()
        
    def get_spacegroup(self, symprec=None, symbol_type=None):
        if symprec == None: symprec = self.symprec = 1e-5
        return spglib.get_spacegroup(self.cell, symprec=symprec, symbol_type=symbol_type) 

    def get_symmetry(self, cell=None, symprec=None): 
        if cell == None: cell = self.cell
        if symprec == None: symprec = self.symprec = 1e-5
        sym = spglib.get_symmetry(cell, symprec)
        if sym != None:
            return sym['rotations'], sym['translations'], sym['equivalent_atoms']  
        else:
            return None
        
    def refine_cell(self, cell=None, symprec=None):     
        if cell == None: cell = self.cell
        if symprec == None: symprec = self.symprec = 1e-5
        cell = spglib.refine_cell(cell, symprec)
        if cell != None:
            return cell
        else:
            return None
        
    def find_primitive(self, cell=None, symprec=None):    
        if cell == None: cell = self.cell
        if symprec == None: symprec = self.symprec = 1e-5
        cell = spglib.find_primitive(cell, symprec)
        if cell != None:
            return cell
        else:
            return None   
  

    def standardize_cell(self, cell=None, to_primitive=False, no_idealize=False, symprec=None):    
        if cell == None: cell = self.cell
        if symprec == None: symprec = self.symprec = 1e-5
        cell = spglib.standardize_cell(cell, to_primitive, no_idealize, symprec)
        if cell != None:
            return cell
        else:
            return None
        
    def get_symmetry_dataset(self, cell=None, symprec=None, angle_tolerance=-1.0, hall_number=0):    
        if cell == None: cell = self.cell
        if symprec == None: symprec = self.symprec = 1e-5
        dataset = spglib.get_symmetry_dataset(cell, symprec, angle_tolerance, hall_number)
        if dataset != None:
            self.number = dataset['number'] 
            self.international = dataset['international'] 
            self.hall = dataset['hall'] 
            self.hall_number = dataset['hall_number'] 
            self.choice = dataset['choice']
            self.transformation_matrix = dataset['transformation_matrix']
            self.origin_shift = dataset['origin_shift']
            self.wyckoffs = dataset['wyckoffs']
            self.site_symmetry_symbols = dataset['site_symmetry_symbols']
            self.equivalent_atoms = dataset['equivalent_atoms']
            self.mapping_to_primitive = dataset['mapping_to_primitive']
            self.rotations = dataset['rotations']    
            self.translations = dataset['translations']   
            self.pointgroup = dataset['pointgroup']  
            self.std_cell = [dataset['std_lattice'], dataset['std_positions'], dataset['std_types']] 
            self.std_rotation_matrix = dataset['std_rotation_matrix']  
            self.std_mapping_to_primitive = dataset['std_mapping_to_primitive']
        else:
            return None

    def get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0]):
        out = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=is_shift)
        if out != None:
            return out[0], out[1]
        else: 
            return None

        
def get_symmetry_from_database(hall_number):    
    dataset = spglib.get_symmetry_from_database(hall_number)
    if dataset != None:
        return dataset['rotations'], dataset['translations']
    else:
        return None
        
def get_spacegroup_type(hall_number):  
    spacegroup = spglib.get_spacegroup_type(hall_number)
    for ele in spacegroup:
        print('%s : %s' % (ele, spacegroup[str(ele)]))
    
def get_hall_number_from_symmetry(rotations, translations, symprec=1e-5):
    return spglib.get_hall_number_from_symmetry(rotations, translations, symprec)
    
def niggli_reduce(lattice, eps=1e-5):
    return spglib.niggli_reduce(lattice, eps)
    
def delaunay_reduce(lattice, eps=1e-5):
    return spglib.delaunay_reduce(lattice, eps)
    
    
        
        

        

        