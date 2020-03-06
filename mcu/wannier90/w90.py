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

import os
import numpy as np
from ..wannier90 import w90_io
from ..utils import plot
from ..cell import spg_wrapper, cell_io
from ..cell import utils as cell_utils
import matplotlib as mpl
import matplotlib.pyplot as plt
        
class main:
    def __init__(self,  seedname="wannier90"):
        '''
            path        : the project directory
            vaspruns    : a str or a list of string as names for *.xml files
            outcars     : a str or a list of string as names for OUTCAR files
        '''
        self.w90_io = w90_io.io(seedname)
        self.w90_io.read_win()
        self.cell = self.w90_io.cell
        self.atom = self.w90_io.atom
        self.kpts = self.w90_io.kpts
        self.band = None

############ Symmetry #################      
    def get_symmetry(self, cell=None, symprec=1e-5, print_atom=False):
        '''Get space group information'''
        if cell == None: 
            cell = self.cell
            is_std, is_prim = spg_wrapper.get_sym(cell, symprec, print_atom)
            self.cell_type = [is_std, is_prim]
        else:
            is_std, is_prim = spg_wrapper.get_sym(cell, symprec)
        
    def to_convcell(self, cell=None, symprec=1e-5):
        '''Transform the unit cell to the standard cell'''
        if cell == None: 
            cell = self.cell
            self.cell = spg_wrapper.cell_to_std(cell, symprec)
            self.cell_type[0] = True
        else:
            return spg_wrapper.cell_to_std(cell, symprec)
            
    def to_primcell(self, cell=None, symprec=1e-5):
        '''Transform the unit cell to the primitive cell'''
        if cell == None: 
            cell = self.cell
            self.cell = spg_wrapper.cell_to_prim(cell, symprec)
            self.cell_type[1] = True
        else:
            return spg_wrapper.cell_to_prim(cell, symprec)      

    def write_poscar(self, cell=None, filename=None):
        if cell == None: cell = self.cell
        cell_io.write_poscar(cell, filename)
        
    def write_cif(self, cell=None, symprec=1e-5, filename=None, symmetry=True):
        if cell == None: 
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
        if cell == None: cell = self.cell
        cell_io.write_xsf(cell, filename) 

############ Plotting ################# 
    def get_efermi(self, num_vb):
        '''E_fermi is assumed to be the valence band maximum. This is a reasonable estimation for insulators'''
        if self.band is None:
            self.w90_io.read_band()
            self.band = self.w90_io.band
        VBM = self.band[:,:num_vb].max()
        return VBM
               
    def get_bandgap(self, efermi=None):
        '''Get the bandgap'''
        assert efermi is not None, "you need to provide the Fermi energy or estimate it using the get_efermi function"
            
        if self.band is None:
            self.w90_io.read_band()
            self.band = self.w90_io.band
        
        CBM = None
        nbands = self.band.shape[1]
        for bandth in range(nbands):
            shifted_band = self.band[:,bandth] - efermi
            if (shifted_band > 0.0).all() == True:
                CBM = self.band[:, bandth]
                VBM = self.band[:, bandth -1]                
                break
            elif ((shifted_band < 0.0).any() == True) and ((shifted_band > 0.0).any() == True):
                print("This is a metal")
                break
                
        if CBM is not None:
            vbm_idx = np.argmax(VBM)
            cbm_idx = np.argmin(CBM)
            bandgap = CBM[cbm_idx] - VBM[vbm_idx]
            direct = False
            if vbm_idx == cbm_idx: direct = True
            
            kpath_frac = self.w90_io.kpath_frac
            print('E(VBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (VBM[vbm_idx], 
                                                            kpath_frac[vbm_idx,0], kpath_frac[vbm_idx,1], kpath_frac[vbm_idx,2]))
            print('E(CBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (CBM[cbm_idx], 
                                                            kpath_frac[cbm_idx,0], kpath_frac[cbm_idx,1], kpath_frac[cbm_idx,2]))
            if direct == True: 
                print('Direct bandgap   : %6.3f' % (bandgap))             
            else:  
                print('Indirect bandgap : %6.3f' % (bandgap))              
                gap1 = CBM[vbm_idx] - VBM[vbm_idx]
                gap2 = CBM[cbm_idx] - VBM[cbm_idx]
                direct_gap = min(gap1, gap2)
                print('Direct bandgap   : %6.3f' % (direct_gap))
 
    def _generate_band(self, efermi=None, spin=0, label=None):
        '''Processing/collecting the band data before the plotting function
           TODO: spin != 0 case will be updated later
        '''
        if self.band is None:
            self.w90_io.read_band()
            self.band = self.w90_io.band
        if efermi is None: efermi = 0.0
            
        assert self.w90_io.klabel is not None, "Cannot find the label for high symmetric k-point in *.win file"
        
        band = self.w90_io.band - efermi
        klabel = self.w90_io.klabel
        label = []
        frac_kpts = [] 
        for kpt in klabel:
            label.append(kpt[0])
            frac_kpts.append(kpt[1])
            
        a = self.cell[0]                        # row vectors
        b = 2*np.pi*np.linalg.inv(a).T     # row vectors
        frac_kpts = np.asarray(frac_kpts)
        abs_kpts = frac_kpts.dot(b)   
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
        
        return band, self.w90_io.proj_kpath, sym_kpoint_coor, label
        
    def plot_band(self, efermi=None, spin=0, save=False, band_color=['#007acc','#808080','#808080'],
                    figsize=(6,6), figname='BAND', xlim=None, ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                                  0 or 1 for spin polarized
                color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                                  
                                  
        '''
        assert isinstance(band_color,list)
        assert len(band_color) == 3
        plot.plot_band(self, efermi=efermi, spin=spin, save=save, band_color=band_color,
                figsize=figsize, figname=figname, xlim=xlim, ylim=ylim, fontsize=fontsize, dpi=dpi, format=format)
        
