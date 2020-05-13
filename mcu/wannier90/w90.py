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
from ..utils import plot
from ..utils.misc import check_exist
from ..cell import utils as cell_utils
from ..cell import cell
from . import w90_io

        
class main(cell.main, plot.main):
    def __init__(self,  seedname="wannier90"):
        '''
            path        : the project directory
            vaspruns    : a str or a list of string as names for *.xml files
            outcars     : a str or a list of string as names for OUTCAR files
        '''
        self.seedname = seedname
        self.get_info()

############ General #################
    def get_info(self, filename=None):    
        '''Extract basis information from the vasprun.xml'''
        if filename is None:
            if check_exist(self.seedname + ".win"):
                filename = self.seedname + ".win"
            else:
                assert 0, "Cannot find " + self.seedname + ".win file"
                
        data = w90_io.read_win(filename)
        lattice = data['unit_cell']
        self.atom = data['atom']
        self.element = list(dict.fromkeys(self.atom ))
        self.kpts = data['kpts']
        self.kpath = data['kpath']
        
        # Make a cell object in the spglib format
        if data['frac_coords'] is not None:
            positions = data['frac_coords']
        elif data['abs_coords'] is not None:
            positions = data['abs_coords'] @ np.linalg.inv(lattice)
        numbers = cell_utils.convert_atomtype(self.atom)
        self.cell = (lattice, positions, numbers)

############ Plotting ################# 
    def get_band(self, filename=None):
        '''Get the band structure info'''
        if filename is None:
            if check_exist(self.seedname + "_band.dat"):
                filename = self.seedname + "_band.dat"
                assert check_exist(self.seedname + "_band.kpt"), "Cannot find any seedname_band.kpt file"
            else:
                assert 0, "Cannot find " + self.seedname + "_band.dat file"
                
        proj_kpath, band = w90_io.read_band(filename)
        nkpts, nbands = band.shape
        band = band.reshape(-1, nkpts, nbands)
        kpath_frac = w90_io.read_kpt(self.seedname + "_band.kpt")
        
        return kpath_frac, proj_kpath, band
                
    def get_efermi(self, num_vb, filename=None):
        '''E_fermi is assumed to be the valence band maximum. This is a reasonable estimation for insulators
           band dimension = [spin, kpt, band]
        '''
        if filename is None:
            if check_exist(self.seedname + "_band.dat"):
                filename = self.seedname + "_band.dat"
            else:
                assert 0, "Cannot find any seedname_band.dat file"
                
        kpath_frac, proj_kpath, band = self.get_band()
        VBM = band[:,:,:num_vb].max()
        return VBM
               
    def get_bandgap(self, efermi=None, spin=0, filename=None):
        '''Get the bandgap'''
        assert efermi is not None, "you need to provide the Fermi energy or estimate it using the get_efermi function"
                
        kpath_frac, proj_kpath, band = self.get_band()
        nspin, nkpts, nbands = band.shape
        for spin in range(nspin):
            print('Spin:', spin)  
            CBM = None
            for bandth in range(nbands):
                shifted_band = band[spin,:,bandth] - efermi
                if (shifted_band > 0.0).all() == True:
                    CBM = band[spin,:,bandth]
                    VBM = band[spin,:,bandth -1]                
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
 
    def _generate_band(self, efermi=None, spin=0, klabel=None):
        '''Processing/collecting the band data before the plotting function
           TODO: spin != 0 case will be updated later
        '''
        if efermi is None: efermi = 0.0
        kpath_frac, proj_kpath, band = self.get_band()
        band = band[spin] - efermi
        
        # Process the labels for k-point
        assert self.kpath is not None, "Cannot find the label for high symmetric k-point in *.win file"
        k1 = self.kpath[0][0]
        k2 = self.kpath[0][1]
        klabel = [k1[0], k2[0]]
        frac_kpts = [k1[1], k2[1]] 
        for i, path in enumerate(self.kpath[1:]):
            k1 = path[0]
            k2 = path[1]
            if (k1[0] != klabel[i + 1]):
                klabel[i + 1] = klabel[i + 1] + "|" + k1[0]
            klabel.append(k2[0])
            frac_kpts.append(k2[1])

        a = self.cell[0]                        # row vectors
        b = 2*np.pi*np.linalg.inv(a).T     # row vectors
        abs_kpts = np.asarray(frac_kpts).dot(b)   
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
        
        return band, proj_kpath, sym_kpoint_coor, klabel
        