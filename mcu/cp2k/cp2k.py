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
from ..utils import plot
from ..utils.misc import check_exist
from ..cell import spg_wrapper, cell_io, cell
from ..cell import utils as cell_utils
from . import cp2k_io
        
class main(cell.main, plot.main):
    def __init__(self,  prefix=None):
        '''
            Initiate stuff, needed to write ...
        '''
        assert prefix is not None, "Provide a prefix name for your project, for example, prefix.out, prefix.bs, etc."
        self.prefix = prefix
        self.get_info() 

############ General #################
        
    def get_info(self, filename=None):    
        '''Extract basis information from the vasprun.xml'''
        if filename is None:
            if check_exist(self.prefix + ".out"):
                filename = self.prefix + ".out"
            else:
                assert 0, "Cannot find any prefix.out file"
                
        data = cp2k_io.read_out(filename)
        self.cell = data['cell']
        self.atom = data['atom']
        self.nelec = int(sum(data['Zeff']))
        
############ Plotting ################# 
    def get_efermi(self, filename=None, set_block=-1):
        '''E_fermi is assumed to be the valence band maximum. This is a reasonable estimation for insulators
        Note: CP2K keeps writing band structure information to the same *.bs file.
        Hence, *.bs file can contain many blocks, each block may have different number of bands and kpts. By default, the last block is used but user can change it by set set_block
        '''
        if filename is None:
            if check_exist(self.prefix + ".bs"):
                filename = self.prefix + ".bs"
            else:
                assert 0, "Cannot find any prefix.bs file"

        data = cp2k_io.read_bs(filename)
        band = data['band'][set_block]
        nspin, nkpts, nbands = band.shape
        
        # Determine efermi
        nelectron = nkpts * self.nelec
        nocc = nelectron // 2       #TODO: uhf is considered 
        energy = band.flatten()
        idx = energy.argsort()
        efermi = energy[idx][:nocc].max()
        
        return efermi

    def get_bandgap(self, filename=None, set_block=-1, efermi=None):
        '''Get the bandgap
        Note: CP2K keeps writing band structure information to the same *.bs file.
        Hence, *.bs file can contain many blocks, each block may have different number of bands and kpts. By default, the last block is used but user can change it by set set_block
        ''' 
        if filename is None:
            if check_exist(self.prefix + ".bs"):
                filename = self.prefix + ".bs"
            else:
                assert 0, "Cannot find any prefix.bs file"
                
        if efermi is None: 
            efermi = self.get_efermi(filename, set_block)
            
        data = cp2k_io.read_bs(filename)
        band = data['band'][set_block]
        kpath_frac = data['kpath_frac'][set_block]
        
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
 
    def _generate_band(self, filename=None, set_block=-1, efermi=0.0, spin=0, klabel=None):
        '''Processing/collecting the band data before the plotting function
           TODO: spin != 0 case will be updated later
        Note: CP2K keeps writing band structure information to the same *.bs file.
        Hence, *.bs file can contain many blocks, each block may have different number of bands and kpts. By default, the last block is used but user can change it by set set_block
        '''

        if filename is None:
            if check_exist(self.prefix + ".bs"):
                filename = self.prefix + ".bs"
            else:
                assert 0, "Cannot find any prefix.bs file"
            
        if efermi is None: 
            efermi = self.get_efermi(filename, set_block)
            
        data = cp2k_io.read_bs(filename)
        band = data['band'][set_block]
        kpath_frac = data['kpath_frac'][set_block]
        symm_k_coor = data['symm_k_coor'][set_block]
        klabel = data['symm_k_label'][set_block]        
        
        # Find absolute kpts
        lattice = self.cell[0]
        b =  2*np.pi*np.linalg.inv(lattice).T               # Get the reciprocal lattice
        abs_kpts = kpath_frac.dot(b)                  # From fractional to absolute in A^-1 unit
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        path = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            
        band = band[spin] - efermi
        a = self.cell[0]                        # row vectors
        b = 2*np.pi*np.linalg.inv(a).T     # row vectors
        abs_kpts = symm_k_coor.dot(b)   
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
        
        return band, path, sym_kpoint_coor, klabel
