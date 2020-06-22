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
from . import cell_io, spg_wrapper
from . import utils as cell_utils
from ..utils import misc


class main:
    def __init__(self, cell=None, cif=None):
        if cell is not None:
            self.cell = cell
        if cif is not None:
            self.cell = self.cif2cell(cif)
         
    def get_symmetry(self, cell=None, symprec=1e-3, verbose='short', angle_tolerance=-1.0, hall_number=0):
        '''Get space group information'''
        if cell is None: cell = self.cell
        spg_wrapper.get_sym(cell, symprec, verbose, angle_tolerance, hall_number)
        
    def to_std_cell(self, cell=None, symprec=1e-3):
        '''Transform the unit cell to the standard cell'''
        if cell is None: cell = self.cell
        return spg_wrapper.cell_to_std(cell, symprec)
            
    def to_prim_cell(self, cell=None, symprec=1e-3):
        '''Transform the unit cell to the primitive cell'''
        if cell is None: cell = self.cell
        return spg_wrapper.cell_to_primitive(cell, symprec)      

    def cif2cell(self, filename):
        '''Return a cell object from a CIF'''
        return cell_io.cif2cell(filename) 
        
    def cif2poscar(self, filename):
        '''Return a cell object from a CIF'''
        cell = cell_io.cif2cell(filename) 
        return self.write_poscar(cell)
        
    def cif2xsf(self, filename):
        '''Return a cell object from a CIF'''
        cell = cell_io.cif2cell(filename) 
        return self.write_xsf(cell)
        
    def cif2pyscf(self, filename, cell_type='primitive'):
        '''Return a PySCF cell object from a CIF
        Input:
        =====
             - CIF
             - gto module from PySCF  
            
        '''
        cell = cell_io.cif2cell(filename) 
        if cell_type == 'primitive':
            lattice, frac_coords, atom_Z = self.to_prim_cell(cell) 
        elif cell_type == 'standard':
            lattice, frac_coords, atom_Z = self.to_std_cell(cell)
        else:
            lattice, frac_coords, atom_Z = cell

        atom_symbol = cell_utils.convert_atomtype(atom_Z)
        abs_coods = frac_coords @ lattice
        atoms = [[atm, abs_coods[i]] for i, atm in enumerate(atom_symbol)]
        return lattice, atoms

    def write_poscar(self, cell=None, filename=None):
        '''Export POSCAR from a cell object'''
        if cell is None: cell = self.cell
        if filename is None: filename = 'POSCAR_by_mcu'
        cell_io.write_poscar(cell, filename)
        
    def write_xsf(self, cell=None, filename=None):
        '''Export xsf from a cell object'''
        if cell is None: cell = self.cell
        if filename is None: filename = 'xsf_by_mcu'
        cell_io.write_xsf(cell, filename) 
        
    def write_cif(self, cell=None, symprec=1e-3, filename=None, symmetrize=True, angle_tolerance=-1.0, hall_number=0):
        '''Export CIF from a cell object'''
        if cell is None: cell = self.cell
        if filename is None: filename = 'cif_by_mcu'
        
        if symmetrize is True:  
            spacegroup, irred_cell, rotations, translations = \
            spg_wrapper.get_sym(cell, symprec, verbose=None, angle_tolerance=angle_tolerance, hall_number=hall_number)
            misc.print_msg("A symmetrized structure (No. {:d}, {:s}) is written in {:s}.cif".format(spacegroup[0], spacegroup[1], filename))
        else:
            spacegroup = [1,'P1']
            irred_cell = cell
            symopt = spg_wrapper.get_symmetry_from_database(1)
            rotations = symopt['rotations']
            translations = symopt['translations']
            misc.print_msg("A P1 structure is written in {:s}.cif".format(filename))
            
        symopt = cell_utils.symop_mat2xyz(rotations, translations)
        cell_io.write_cif(irred_cell, spacegroup, symopt, filename) 
        
    def get_mapping_kpts(self, cell=None, mesh=[2,2,2], is_shift=[0, 0, 0], symprec=1e-3, is_time_reversal=True, no_spatial=False):
        '''Generate a uniform k-mesh taking into account the space group and time-reversal symmetry'''
        if cell is None: cell = self.cell
        if no_spatial:
            cell = (np.random.rand(3,3), cell[1], cell[2])
            mapping, grid = spg_wrapper.get_ir_reciprocal_mesh(mesh=mesh, cell=cell, is_shift=is_shift, \
                        symprec=symprec, is_time_reversal=is_time_reversal)
        else:
            mapping, grid = spg_wrapper.get_ir_reciprocal_mesh(mesh=mesh, cell=cell, is_shift=is_shift, \
                        symprec=symprec, is_time_reversal=is_time_reversal)
                        
        return mapping, grid
        
