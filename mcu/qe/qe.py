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
from ..vasp import const
from . import qe_io

        
class main:
    def __init__(self,  prefix=None):
        '''
        
        '''
        if prefix is None:
            print("Provide a prefix name for your project, for example, prefix.scf.out, prefix.band.out, etc. Otherwise you would have to specify each output file when performing analysis")
        self.prefix = prefix
        

############ Plotting ################# 
    def get_efermi(self, data=None):
        '''Get Fermi energy'''
        if data is None:
            if check_exist(self.prefix + ".band.out"):
                filename = self.prefix + ".band.out"
            elif check_exist(self.prefix + ".scf.out"):
                filename = self.prefix + ".scf.out"
            data = qe_io.read_pw_output(filename)
            
        if data['hoco'] is not None:
            efermi = data['hoco']
        elif data['efermi'] is not None:
            efermi = data['efermi']
        else:
            # Take the HOMO as the E Fermi
            band = data['eigenvals']
            nelec = int(data['nelec']) 
            if band.shape[0] == 1:
                npairs = nelec//2
                efermi = band[0,:,:npairs].max()
            else:  
                nspin, nkpts, nband = band.shape
                band = np.sort(band.flatten())
                nelecs = nelec * nkpts
                efermi = band[:nelecs].max()
        
        return efermi


    def get_band(self, filename=None):
        '''
        if filename is not specified, band structure is obtained by searching for:
            (1) self.prefix + ".band.out"
            (2) self.prefix + ".scf.out"
        '''
        if filename is None:
            if check_exist(self.prefix + ".band.out"):
                filename = self.prefix + ".band.out"
            elif check_exist(self.prefix + ".scf.out"):
                filename = self.prefix + ".scf.out"
            else:
                assert 0, "Cannot find any band structure file"
                
        # Collecting data from pw_output
        data = qe_io.read_pw_output(filename)
        band = data['eigenvals']
        kpts = data['kpts'][:,:3]       # the last column is the kpt weight
        alat = data['alat'] * const.AUTOA       # in angstrom
        lattice = data['crystal axes'] * alat
        recip_lattice = data['reciprocal axes'] * 2*np.pi / alat
        
        # Find absolute kpts
        abs_kpts = kpts.dot(recip_lattice)                  # From fractional to absolute in A^-1 unit
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
        
        efermi = self.get_efermi(data)
        
        return band, kpts, proj_kpath, recip_lattice, efermi

    def get_bandgap(self, filename=None):
        '''Get the bandgap'''
        
        band, kpath_frac, proj_kpath, recip_lattice, efermi = self.get_band(filename) 
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
 
    def _generate_band(self, filename=None, efermi=None, spin=0, label=None):
        '''Processing/collecting the band data before the plotting function
           label            : a list of labels and corresponding coordinates
        '''
        band, kpath_frac, proj_kpath, recip_lattice, efermi_ = self.get_band(filename)
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
        # Find absolute coordinates for high symmetric kpoints  
        sym_kpoint_coor = None
        if label is not None:
            assert isinstance(label,list)           # label needs to be a list of labels and corresponding coordinates
            temp = []
            coor_kpts = [] 
            for kpt in label:
                temp.append(kpt[0])
                coor_kpts.append(kpt[1:])
            label = temp       
            coor_kpts = np.asarray(coor_kpts)
            abs_kpts = coor_kpts.dot(recip_lattice)   
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
                
        return band, proj_kpath, sym_kpoint_coor, label
        
    def plot_band(self, efermi=None, label=None, spin=0, save=False, band_color=['#007acc','#808080','#808080'],
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
                figsize=figsize, figname=figname, xlim=xlim, ylim=ylim, fontsize=fontsize, dpi=dpi, format=format, label=label)
                

        
