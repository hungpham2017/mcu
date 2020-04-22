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
import matplotlib as mpl
import matplotlib.pyplot as plt

      
      
def plot_band(self, efermi=None, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'],
                figsize=(6,6), figname='BAND', xlim=None, ylim=[-6,6], fontsize=18, dpi=300, format='png'):
    '''Plot band structure
       
        Attribute:
            efermi          : a Fermi level or a list of Fermi levels
            spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                              0 or 1 for spin polarized
            color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                              
                              
    '''
    
    assert isinstance(band_color,list)
    assert len(band_color) == 3
    
    band, path, sym_kpoint_coor, label = self._generate_band(efermi=efermi, spin=spin, label=label)  

    ##----------------------------------------------------------
    ##Plotting:        
    ##----------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    yrange = (-500,500)
               
    # Plot the high symmetric kpoint grid
    
    if sym_kpoint_coor is not None:
        for kpt in range(sym_kpoint_coor.shape[0]):
            ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color=band_color[1],linewidth=1.0)
                
    if label is not None and xlim == None:
        nkpts = len(label)
        assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
        for kpt in range(nkpts):   
            point = label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(sym_kpoint_coor[kpt]/path.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)    

    # Plot bands            
    ax.plot([0,path.max()],[0,0],color=band_color[2],linewidth=1.0, dashes=[6,3])       # Fermi level
    for ith in range(band.shape[1]):
        ax.plot(path.T,band[:,ith],color=band_color[0],linewidth=1.0)    
         
    # Graph adjustments             
    ax.tick_params(labelsize=fontsize)
    if xlim == None:
        plt.xlim([0,path.max()])
        plt.xticks([])
        plt.xlabel('k', size=fontsize+4)
    else:
        plt.xlim(xlim)
        plt.xlabel('k ' + r'($\AA^{-1}$)', size=fontsize+4)
    ax.xaxis.set_label_coords(0.5, -0.08) 
    plt.ylabel('Energy (eV)', size=fontsize+4)        
    plt.ylim(ylim)
    plt.tight_layout()
    if save == True: 
        fig.savefig(figname+'.'+format,dpi=dpi,format=format)      
    else:
        plt.show()  
        
def plot_phononband(self, unit='CM', gamma_correct=False, threshold=8.06554, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'],
                figsize=(6,6), figname='PHONONBAND', xlim=None, ylim=None, fontsize=18, dpi=300, format='png'):
    '''Plot phnon band structure
       
        Attribute:
            efermi          : a Fermi level or a list of Fermi levels
            spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                              0 or 1 for spin polarized
            color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                              
                              
    '''
    
    if ylim is None:
        if unit.lower() == "cm": ylim = [-100, 600]
        if unit.lower() == "thz": ylim = [-5, 20]
        if unit.lower() == "mev": ylim = [-15, 80]        
        
    assert isinstance(band_color,list)
    assert len(band_color) == 3
    
    band, path, sym_kpoint_coor, label = self._generate_phononband(unit=unit, gamma_correct=gamma_correct, threshold=threshold, spin=spin, label=label) 

    ##----------------------------------------------------------
    ##Plotting:        
    ##----------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    yrange = (-10000,10000)
               
    # Plot the high symmetric kpoint grid
    for kpt in range(sym_kpoint_coor.shape[0]):
        ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color=band_color[1],linewidth=1.0)
                
    if label is not None and xlim == None:
        nkpts = len(label)
        assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
        for kpt in range(nkpts):   
            point = label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(sym_kpoint_coor[kpt]/path.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)    

    # Plot bands            
    ax.plot([0,path.max()],[0,0],color=band_color[2],linewidth=1.0, dashes=[6,3])       # Fermi level
    for ith in range(band.shape[1]):
        ax.plot(path.T,band[:,ith],color=band_color[0],linewidth=1.0)    
         
    # Graph adjustments             
    ax.tick_params(labelsize=fontsize)
    if xlim == None:
        plt.xlim([0,path.max()])
        plt.xticks([])
        plt.xlabel('Wave vector', size=fontsize+2)
    else:
        plt.xlim(xlim)
        plt.xlabel('Wave vector', size=fontsize+2)
    ax.xaxis.set_label_coords(0.5, -0.08) 
    if unit.lower() == 'thz':
        plt.ylabel('Frequency (THz)', size=fontsize+2)  
    elif unit.lower() == 'cm':
        plt.ylabel('Frequency ' + r'(cm$^{-1}$)', size=fontsize+2)  
    elif unit.lower() == 'mev':
        plt.ylabel('Frequency (meV)', size=fontsize+2)     

    plt.ylim(ylim)
    plt.tight_layout()
    if save == True: 
        fig.savefig(figname+'.'+format,dpi=dpi,format=format)      
    else:
        plt.show() 