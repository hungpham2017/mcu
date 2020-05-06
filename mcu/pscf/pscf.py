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
from ..utils import plot, str_format
from ..utils.misc import check_exist
from ..cell import utils as cell_utils
from ..cell import cell
from ..vasp import const

        
class main(cell.main, plot.main):
    def __init__(self,  pyscf_cell):
        self.pyscf_cell = pyscf_cell
        self.get_info()
        
############ General #################
        
    def get_info(self):    
        '''Extract basis information from the vasprun.xml'''
                
        self.nelec = self.pyscf_cell.nelectron
        self.nao = self.pyscf_cell.nao
        
        # Make a cell object in the spglib format
        atom = []
        abs_positions = []
        for atm in self.pyscf_cell.atom:
            atom.append(atm[0])
            abs_positions.append(atm[1])
            
        self.atom  = atom
        self.natom  = len(self.atom)
        self.element = list(dict.fromkeys(self.atom))
        lattice = self.pyscf_cell.lattice_vectors()*const.AUTOA
        positions = abs_positions @ np.linalg.inv(lattice)
        numbers = cell_utils.convert_atomtype(self.atom)
        self.cell = (lattice, positions, numbers)
        self.band = None
        self.kpts = None        
        
############ Plotting ################# 
    def make_kpts(self, kpath, npoint=20):
        '''Giving a 1D k-path in fractional coordinates, return the absolute kpts list'''
        labels, coords = str_format.format_klabel(kpath)
        frac_kpts = []
        for path_th in range(len(coords) - 1):
            temp = np.arange(npoint)/(npoint - 1)
            path = np.ones([npoint,3]) * np.asarray([temp]*3).T * (coords[path_th + 1] - coords[path_th]) + coords[path_th]    
            if path_th != len(coords) - 2:
                path = path[:-1,:]
            frac_kpts.append(path)

        frac_kpts = np.concatenate(frac_kpts)
        return frac_kpts       
        
    def get_band(self, band=None, kpts=None):
        if band is None: 
            assert self.band is not None, "You need to set mcu.bands"
            band = np.asarray(self.band[0])
            nkpts = band.shape[-2]
            nband = band.shape[-1]
            band = band.reshape(-1, nkpts, nband) * const.AUTOEV
        if kpts is None: 
            assert self.kpts is not None, "You need to set mcu.kpts"
            kpts = self.kpts     
            
        kpath_frac = kpts
        
        # Find absolute kpts and shift the band
        lattice = self.cell[0]
        recip_lattice = 2 * np.pi * np.linalg.inv(lattice).T # Get the reciprocal lattice in the row vector format
        abs_kpts = kpts.dot(recip_lattice)                  # From fractional to absolute
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            
        # Determine efermi
        nelectron = nkpts * self.nelec
        nocc = nelectron // 2       #TODO: uhf is considered 
        energy = band.flatten()
        idx = energy.argsort()
        efermi = energy[idx][:nocc].max()
        
        return band, kpath_frac, proj_kpath, recip_lattice, efermi
    
    def _generate_band(self, efermi=None, spin=0, klabel=None):
        '''Processing/collecting the band data before the plotting function
           klabel            : a list of labels and corresponding coordinates for high symmetry k-points
        '''          
        band, kpath_frac, proj_kpath, recip_lattice, efermi_ = self.get_band()
        
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
        # Find absolute coordinates for high symmetric kpoints  
        sym_kpoint_coor = None
        if klabel is not None:
            klabel, coor_kpts = str_format.format_klabel(klabel)  
            assert coor_kpts is not None, "You need to provide the coordinates for high symmetric k in the klabel"  
            abs_kpts = coor_kpts.dot(recip_lattice)   
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
                
        return band, proj_kpath, sym_kpoint_coor, klabel

    def get_bandgap(self):
        '''Get the bandgap'''
                
        band, kpath_frac, proj_kpath, recip_lattice, efermi = self.get_band()
        nspin, nkpts, nbands = band.shape
        for spin in range(nspin):
            print('Spin:', spin)  
            CBM = None
            for bandth in range(nbands):
                shifted_band = band[spin,:,bandth] - efermi
                if (shifted_band > 0.0).all() == True:
                    CBM = band[spin,:, bandth]
                    VBM = band[spin,:, bandth -1]                
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

                print('  E(VBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (VBM[vbm_idx], 
                                                                kpath_frac[vbm_idx,0], kpath_frac[vbm_idx,1], kpath_frac[vbm_idx,2]))
                print('  E(CBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (CBM[cbm_idx], 
                                                                kpath_frac[cbm_idx,0], kpath_frac[cbm_idx,1], kpath_frac[cbm_idx,2]))
                if direct == True: 
                    print('  Direct bandgap   : %6.3f' % (bandgap))             
                else:  
                    print('  Indirect bandgap : %6.3f' % (bandgap))              
                    gap1 = CBM[vbm_idx] - VBM[vbm_idx]
                    gap2 = CBM[cbm_idx] - VBM[cbm_idx]
                    direct_gap = min(gap1, gap2)
                    print('  Direct bandgap   : %6.3f' % (direct_gap))
 
