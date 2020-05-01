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
from . import cell_io, spg_wrapper
        
        
class main:
    def __init__(self, cell):
        self.cell = cell
         
    def get_symmetry(self, cell=None, symprec=1e-5, print_atom=False):
        '''Get space group information'''
        if cell is None: 
            cell = self.cell
            print("HOHOH", spg_wrapper.get_sym(cell, symprec, print_atom))
            is_std, is_prim = spg_wrapper.get_sym(cell, symprec, print_atom)
            self.cell_type = [is_std, is_prim]
        else:
            is_std, is_prim = spg_wrapper.get_sym(cell, symprec)
        
    def to_convcell(self, cell=None, symprec=1e-5):
        '''Transform the unit cell to the standard cell'''
        if cell is None: 
            cell = self.cell
            self.cell = spg_wrapper.cell_to_std(cell, symprec)
            self.cell_type[0] = True
        else:
            return spg_wrapper.cell_to_std(cell, symprec)
            
    def to_primcell(self, cell=None, symprec=1e-5):
        '''Transform the unit cell to the primitive cell'''
        if cell is None: 
            cell = self.cell
            self.cell = spg_wrapper.cell_to_prim(cell, symprec)
            self.cell_type[1] = True
        else:
            return spg_wrapper.cell_to_prim(cell, symprec)      

    def write_poscar(self, cell=None, filename=None):
        if cell is None: cell = self.cell
        cell_io.write_poscar(cell, filename)
        
    def write_cif(self, cell=None, symprec=1e-5, filename=None, symmetry=True):
        if cell is None: 
            cell = self.cell
            is_std, is_prim = self.cell_type 
            if is_std and symmetry==True: 
                cell = self.to_stdcell(cell, symprec) 
                spacegroup, equi_atoms, rotations, translations = spg_wrapper.get_sym(cell, symprec, export_operator=True)
            elif is_prim and symmetry==True:
                cell = self.to_primcell(cell, symprec)
                spacegroup, equi_atoms, rotations, translations = spg_wrapper.get_sym(cell, symprec, export_operator=True)
            else:
                spacegroup = ['1','P1']
                equi_atoms = np.arange(len(cell[2]))
                symopt = spg_wrapper.get_symmetry_from_database(1)
                rotations, translations = symopt['rotations'], symopt['translations']
        else:
            spacegroup = ['1','P1']
            equi_atoms = np.arange(len(cell[2]))
            symopt = spg_wrapper.get_symmetry_from_database(1)
            rotations, translations = symopt['rotations'], symopt['translations']
        symopt = cell_utils.symop_mat2xyz(rotations, translations)
        cell_io.write_cif(cell, spacegroup, equi_atoms, symopt, filename) 

    def write_xsf(self, cell=None, filename=None):
        if cell is None: cell = self.cell
        cell_io.write_xsf(cell, filename) 
        

        

        