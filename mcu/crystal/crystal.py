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
from ..vasp import const
from . import crystal_io
import matplotlib as mpl
import matplotlib.pyplot as plt
        
class main:
    def __init__(self,  seedname="outfile"):
        '''
        
        '''
        self.seedname = seedname
        

############ Plotting ################# 
    def get_band(self, phonon=False, gamma_correct=False, threshold=8.06554):
        '''Make a band from from 
           NOTE: the proj_kpath is computed from the dk and nkp. Due to the round out error in f25, the computed high symmetric k-point coordinates won't be exactly the same as values obtained from *.BAND file.
           
           if phonon == True: the phonon band is read. Three small negative acoustic modes at Gamma can be removed by using the gamma_correct keyword and threshold 1 meV = 8.06554 cm^-1  
        '''
        data = crystal_io.read_f25(self.seedname + ".f25")
        temp = []
        ihferm_list = []
        proj_kpath = []
        sym_kpoint_coor = [0.0]
        shift = 0
        for block in data:
            ihferm, type, nband, nkp, dk, efermi, eigenvals = block
            if type == 'BAND':
                temp.append(eigenvals.reshape(-1, nband))
                path = np.arange(nkp)*dk + shift
                proj_kpath.append(path)
                shift = path[-1] + dk
                sym_kpoint_coor.append(path[-1])
        
        if ihferm == 0 or ihferm == 2:
            band = np.float64([np.vstack(temp)])
            proj_kpath = np.hstack(proj_kpath)
        elif ihferm == 1 or ihferm == 3:
            nblock = len(temp) // 2
            band_up = np.vstack(temp[:nblock])            
            band_down = np.vstack(temp[nblock:])  
            band = np.float64([band_up, band_down])
            proj_kpath = np.hstack(proj_kpath[:nblock])
        sym_kpoint_coor = np.float64(sym_kpoint_coor)
        if phonon == True:
            if gamma_correct == True:
                # Should be checked first before using the correction to make sure it is really acoustic phonon modes
                nspin, nkpts, nband = band.shape
                band = band.flatten()
                imag_mode_idx = band < 0.0
                imag_mode = band[imag_mode_idx]
                imag_mode_idx2 = imag_mode > -threshold
                imag_mode[imag_mode_idx2] = 0.0
                band[imag_mode_idx] = imag_mode
                band = band.reshape(nspin, nkpts, nband)
        else:
            band = const.AUTOEV * band
            efermi = const.AUTOEV * efermi
            
        return band, proj_kpath, sym_kpoint_coor, efermi

    def get_bandgap(self):
        '''Get the bandgap'''
        
        band, proj_kpath, sym_kpoint_coor, efermi = self.get_band() 
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
                
                # TODO: the kpath_frac currently cannot be obtained form f25.
                # Other outputs are needed
                # kpath_frac = self.cp2k_io.kpath_frac[set_block]
                # print('  E(VBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (VBM[vbm_idx], 
                                                                # kpath_frac[vbm_idx,0], kpath_frac[vbm_idx,1], kpath_frac[vbm_idx,2]))
                # print('  E(CBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (CBM[cbm_idx], 
                                                                # kpath_frac[cbm_idx,0], kpath_frac[cbm_idx,1], kpath_frac[cbm_idx,2]))
                if direct == True: 
                    print('  Direct bandgap   : %6.3f' % (bandgap))             
                else:  
                    print('  Indirect bandgap : %6.3f' % (bandgap))              
                    gap1 = CBM[vbm_idx] - VBM[vbm_idx]
                    gap2 = CBM[cbm_idx] - VBM[cbm_idx]
                    direct_gap = min(gap1, gap2)
                    print('  Direct bandgap   : %6.3f' % (direct_gap))
 
    def _generate_band(self, efermi=None, spin=0, label=None):
        '''Processing/collecting the band data before the plotting function
        '''
        band, proj_kpath, sym_kpoint_coor, efermi_ = self.get_band()
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
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
                
    def _generate_phononband(self, unit="CM", gamma_correct=False, threshold=8.06554, spin=0, label=None):
        '''Processing/collecting the phnon band (cm^-1) data before the plotting function
           Unit: CM = CM^-1, THZ = THz, MEV = meV
        '''
        band, proj_kpath, sym_kpoint_coor, efermi_ = self.get_band(phonon=True, gamma_correct=gamma_correct, threshold=threshold)        
        if unit.lower() == "thz":
            band = const.CMTOTHZ * band
        elif unit.lower() == "mev":
            band = const.CMTOMEV * band
            
        return band[spin], proj_kpath, sym_kpoint_coor, label
        
    def plot_phononband(self, unit="CM", gamma_correct=False, threshold=8.06554, label=None, spin=0, save=False, band_color=['#007acc','#808080','#808080'],
                    figsize=(6,6), figname='PHONONBAND', xlim=None, ylim=None, fontsize=18, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                                  0 or 1 for spin polarized
                color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                                  
                                  
        '''
        assert isinstance(band_color,list)
        assert len(band_color) == 3
        plot.plot_phononband(self, unit=unit, gamma_correct=gamma_correct, threshold=threshold, spin=spin, save=save, band_color=band_color,
                figsize=figsize, figname=figname, xlim=xlim, ylim=ylim, fontsize=fontsize, dpi=dpi, format=format, label=label)
        
