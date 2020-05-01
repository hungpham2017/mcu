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
from . import str_format
      
      
def plot_band(calculator, efermi=None, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'],
                figsize=(6,6), figname='BAND', xlim=None, ylim=[-6,6], fontsize=18, dpi=300, format='png'):
    '''Plot band structure
       
        Attribute:
            calculator      : the engine to generate band structure (VASP, QE, CRYSTAL, etc.)
            efermi          : a Fermi level or a list of Fermi levels
            spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                              0 or 1 for spin polarized
            color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                              
                              
    '''
    
    assert isinstance(band_color,list)
    assert len(band_color) == 3
    
    band, kpath, sym_kpoint_coor, formatted_label, conventional = calculator._generate_band(efermi=efermi, spin=spin, label=label)  

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
                
    if formatted_label is not None and xlim is None:
        nkpts = len(formatted_label)
        assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
        for kpt in range(nkpts):   
            point = formatted_label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(sym_kpoint_coor[kpt]/kpath.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)    

    # Plot bands            
    ax.plot([0,kpath.max()],[0,0],color=band_color[2],linewidth=1.0, dashes=[6,3])       # Fermi level
    for ith in range(band.shape[1]):
        ax.plot(kpath.T,band[:,ith],color=band_color[0],linewidth=1.0)    
         
    # Graph adjustments             
    ax.tick_params(labelsize=fontsize)
    if xlim is None:
        plt.xlim([0,kpath.max()])
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
        
        
def plot_phononband(calculator, unit='CM', gamma_correct=False, threshold=8.06554, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'], figsize=(6,6), figname='PHONONBAND', xlim=None, ylim=None, fontsize=18, dpi=300, format='png'):
    '''Plot phnon band structure
       
        Attribute:
            calculator      : the engine to generate band structure (VASP, QE, CRYSTAL, etc.)
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
    
    band, kpath, sym_kpoint_coor, formatted_label = calculator._generate_phononband(unit=unit, gamma_correct=gamma_correct, threshold=threshold, spin=spin, label=label) 

    ##----------------------------------------------------------
    ##Plotting:        
    ##----------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    yrange = (-10000,10000)
               
    # Plot the high symmetric kpoint grid
    for kpt in range(sym_kpoint_coor.shape[0]):
        ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color=band_color[1],linewidth=1.0)
                
    if formatted_label is not None and xlim is None:
        nkpts = len(formatted_label)
        assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
        for kpt in range(nkpts):   
            point = formatted_label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(sym_kpoint_coor[kpt]/kpath.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)    

    # Plot bands            
    ax.plot([0,kpath.max()],[0,0],color=band_color[2],linewidth=1.0, dashes=[6,3])       # Fermi level
    for ith in range(band.shape[1]):
        ax.plot(kpath.T,band[:,ith],color=band_color[0],linewidth=1.0)    
         
    # Graph adjustments             
    ax.tick_params(labelsize=fontsize)
    if xlim is None:
        plt.xlim([0,kpath.max()])
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
        
        
def plot_pband(calculator, efermi=None, spin=0, label=None, gradient=False, lm=None, band=None, color=None, band_color=['#007acc','#808080','#808080'], scale=1.0, alpha=0.5, cmap='bwr', edgecolor='none', facecolor=None, marker=None, legend=None, loc="upper right", legend_size=1.0, save=False, figname='pBAND', figsize=(6,6), xlim=None, ylim=[-6,6], fontsize=18, dpi=600, format='png'):
    '''Plot projected band structure
       
        Attribute:
            efermi      : a Fermi level or a list of Fermi levels, it is automatically extracted frim vasprun.xml or OUTCAR
            spin        : 0  for spin unpolarized and LSORBIT = .TRUE.
                              0 or 1 for spin polarized
            label       : label for high symmetric points, e.g. 'G-X-L-W-G',
                              if hybridXC=True, the lavel should be a list of labels plus their coordinates
                              
            ########################PLOTTING STYLE###################################
            Examples for lm:
                lm = 'Ni:s ; p ; d'             :   three groups: (1) s of Ni ; (2) all p orbitals ; (3) all d orbitals
                lm = ['Ni:s,pd', 'O1:p;O2']     :   two groups: (1) s,p,d of Ni ; (2) all p orbitals of the 1st O  and all otbitals of O2 
                lm = ['Ni1;O', 'N']             :   two groups: (1) the 1st Ni and all the O atoms ; (2) All N atom
 
            if gradient == True: user has to provide a TWO groups of orbitals  
                for example, lm = 'Ni:s ; p' or ['Ni:s,pd', 'O1:p;O2']    
            #########################################################################
                         
            band        : the first value indicates the number of valence bands from the VBM
                          the second value indicates the number of conduction bands from the CBM
            color       : a list of strings indicating the color, following matplotlib
            scale       : the size of the marker
            alpha       : the transparency level of curves
            cmap        : color map in the style 3, following the matplotlib
            edgecolor   : the marker's border color in the style 3, default: 'none', any color code should work
            facecolor   : the filling color of style 1 and 2
                          None              : taking from the color arg
                          'none'            : unfilling circle
                          [False,True,...]  : True for filling markers and False for empty ones
            marker      : a list of marker shape, default is: 'o'
            legend      : a list of labels for different group of orbitals (same color) for the style 1 and 2            
    '''
    
    if lm is None:    
        if gradient: 
            lm = 's,p'   
        else:
            lm = [atom+':s,p,d' for atom in calculator.element]  
   
    # Check if the band values are reasonable otherwise generate it
    band_idx = band
    if band_idx is None:
        idx_vbm = int(calculator.nelec)
        if calculator.soc == False: idx_vbm = idx_vbm//2               # Estimation for insulator, generally wrong for metal
        first_band = int(idx_vbm - 5)
        last_band = int(idx_vbm + 4)
        if first_band < 0: first_band = 0
        if last_band > calculator.nbands - 1: last_band = calculator.nbands - 1
        band_idx = [first_band, last_band]
    else:
        assert band_idx[0] <= band_idx[1]                    # from band[0] to band[1]
        if band_idx[0] < 1: band_idx[0] = 1     # index used in OUTCAR, will be shifted to start at zero
        if band_idx[1] > calculator.nbands: band_idx[1] = calculator.nbands              # Cannot larger than the number of bands
        band_idx[0] = band_idx[0] -1
        band_idx[1] = band_idx[1] -1     
    
    band, kpath, sym_kpoint_coor, formatted_label, conventional = calculator._generate_band(efermi=efermi, spin=spin, label=label)  
    pband = calculator._generate_pband(spin=spin, gradient=gradient, lm=lm)
    
    ##----------------------------------------------------------
    ##Plotting:        
    ##----------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    yrange = (-500,500)
    
    # Customization:
    border = 1.08
    
    # Plot the high symmetric kpoint grid
    if sym_kpoint_coor is not None:
        for kpt in range(sym_kpoint_coor.shape[0]):
            ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color=band_color[1],linewidth=1.0)

    if formatted_label is not None and xlim is None:
        nkpts = len(formatted_label)
        assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
        for kpt in range(nkpts):   
            point = formatted_label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(sym_kpoint_coor[kpt]/kpath.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)     
        
    # Plot bands            
    ax.plot([0, kpath.max()], [0,0], color=band_color[2], linewidth=1.0, dashes=[6,3])
    for ith in range(band.shape[1]):
        ax.plot(kpath.T, band[:,ith], color=band_color[0],linewidth=1.0)    
         
    # Plot pbands 
    color_list = ['r','g','b','y','m','c']
    kpath = np.array(kpath).flatten() 
    
    if gradient:
        kpath = np.array(kpath).flatten()
        if marker is None: 
            marker = 'o'
        else:
            assert isinstance(marker,str)
        for ith in range(band_idx[0],band_idx[1]+1):
            plt.scatter(kpath, band[:,ith], c=pband[:,ith], s=50*scale, vmin=0.0, vmax=1., cmap=cmap, marker=marker, edgecolor=edgecolor) 
        cbar = plt.colorbar(ticks=[])
        cbar.outline.set_linewidth(border)
    
    else:
        pband = 200 * scale * np.power(pband,2)     # The radius of the marker ~ the percent 
        
        # Color
        if color is None: 
            color = color_list
        else:
            assert isinstance(color,list) or isinstance(color,str)
            if isinstance(color,str): color = [color]
            
        if facecolor is None: 
            fcolors = color
        elif isinstance(facecolor,list):
            assert len(facecolor) == len(pband)
            fcolors = []
            for i in range(len(facecolor)):
                assert facecolor[i] == True or facecolor[i] == False
                if facecolor[i] == True: fcolors.append(color[i]) 
                if facecolor[i] == False: fcolors.append('none') 
        elif facecolor == 'none':
            fcolors = ['none']*len(pband)

        # Marker
        if marker is None: 
            marker = ['o']*len(pband)
        else:
            assert isinstance(marker,list) or isinstance(legend,str)
            if isinstance(marker,str): marker = [marker]
            assert len(marker) == len(pband)                
        
        # legend    
        if legend is not None:
            assert isinstance(legend,list) or isinstance(legend,str)
            if isinstance(legend,str): legend = [legend]
            assert len(legend) == len(pband)
            
        # Actual plotting
        for lm in range(len(pband)):
            for ith in range(band_idx[0],band_idx[1]):
                ax.scatter(kpath, band[:,ith], s=pband[lm][:,ith], facecolors=fcolors[lm], edgecolors=color[lm], alpha=alpha, marker=marker[lm])
            ith = band_idx[1]
            if legend is None:
                ax.scatter(kpath, band[:,ith], s=pband[lm][:,ith], facecolors=fcolors[lm], edgecolors=color[lm], alpha=alpha, marker=marker[lm])
            else:
                ax.scatter(kpath, band[:,ith], s=pband[lm][:,ith], facecolors=fcolors[lm], edgecolors=color[lm], alpha=alpha, marker=marker[lm],label=legend[lm])                    
            
        if legend is not None: 
            lgnd = ax.legend(loc=loc, numpoints=1, fontsize=fontsize)
            for i in range(len(pband)): lgnd.legendHandles[i]._sizes = [legend_size*60]
            
    
    # Graph adjustments             
    ax.tick_params(labelsize=fontsize, width=border)
    ax.spines['top'].set_linewidth(border)
    ax.spines['right'].set_linewidth(border)
    ax.spines['bottom'].set_linewidth(border)
    ax.spines['left'].set_linewidth(border)
    if xlim is None:
        plt.xlim([0,kpath.max()])
        plt.xticks([])
        plt.xlabel('k', size=fontsize+4)
    else:
        plt.xlim(xlim)
        plt.xlabel('k ' + r'($\AA^{-1}$)', size=fontsize+4)
    ax.xaxis.set_label_coords(0.5, -0.08) 
    if legend is not None: ax.legend()
    plt.ylabel('Energy (eV)', size=fontsize+4)        
    plt.ylim(ylim)
    plt.tight_layout()
    plt.ylabel('Energy (eV)', size=fontsize+4)        
    plt.ylim(ylim)
    plt.tight_layout()
    plt.ylabel('Energy (eV)', size=fontsize+4)        
    plt.ylim(ylim)
    plt.tight_layout()
    if save == True: 
        fig.savefig(figname+'.'+format, dpi=dpi, format=format)      
    else:
        plt.show()
        
        
def plot_dos(calculator, style='horizontal', efermi=None, spin=0, lm=None, color=None, legend=None, loc="upper right", fill=True, alpha=0.2, save=False, figname='DOS', figsize=(6,3), elim=(-6,6), yscale=1.1, fontsize=18, dpi=600, format='png'):
    '''Plot projected band structure
        For multiple vasprun.xml, user can choose one of them to plot the DOS. Default: the first vasprun.xml

        Attribute:
        style           : 'horizontal' or 'vertical'
        
        for VASP only:
            spin            : spin of DOS.
                              For LSORBIT == True: spin = 0,1,2,3
                              For ISPIN = 2      : spin = 0,1
        spin                : updown, visualize both up and down spin
                          
        lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
    '''
        
    
    
    if lm is None: 
        lm = [atom+':s,p,d' for atom in calculator.element] 
        
    formatted_lm = str_format.general_lm(lm)[0]
    if legend is None:
        legend = ['PDOS-' + str(i) for i in range(len(formatted_lm))]
         
    if spin == 'updown':
        if calculator.ispin != 2: raise Exception('ISPIN must be 2 for the up-down DOS plotting')
        tdos0, pdos0 = calculator._generate_dos(efermi=efermi, spin=0, lm=lm)
        tdos1, pdos1 = calculator._generate_dos(efermi=efermi, spin=1, lm=lm)
        tdos1 = -tdos1
        pdos1 = -pdos1
        if figsize == (6,3) and style=='horizontal' : figsize = (6,5)
        if figsize == (6,3) and style=='vertical' : figsize = (5,6)
    else:
        if figsize == (6,3) and style=='vertical' : figsize = (3,6)
        tdos0, pdos0 = calculator._generate_dos(efermi=efermi, spin=spin, lm=lm)
    
    ##----------------------------------------------------------
    ##Plotting:        
    ##----------------------------------------------------------
    color_list = ['k','r','g','b','y','m','c']
    if color is None: color = color_list
    
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    yrange = (-500,500)
    
    # Plot DOS    
    if style == 'horizontal':
        ax.plot(tdos0[:,0], tdos0[:,1], color=color[0],linewidth=1.1,label='TDOS')
        if pdos0 is not None:
            for orb in range(pdos0.shape[1]): 
                ax.plot(tdos0[:,0], pdos0[:,orb], color=color[orb+1], linewidth=1.0, label=legend[orb])
                if fill == True: ax.fill_between(tdos0[:,0], pdos0[:,orb], color=color[orb+1], alpha=alpha)
            
        if spin == 'updown':
            ax.plot(tdos0[:,0], tdos1[:,1], color=color[0], linewidth=1.1)
            if pdos1 is not None:
                for orb in range(pdos1.shape[1]): 
                    ax.plot(tdos0[:,0], pdos1[:,orb], color=color[orb+1], linewidth=1.0)
                    if fill == True: ax.fill_between(tdos0[:,0], pdos1[:,orb], color=color[orb+1], alpha=alpha)
    
        # Graph adjustments 
        plt.xlabel('Energy (eV)', size=fontsize+4)   
        plt.ylabel('DOS', size=fontsize+4)
        if spin == 'updown':
            plt.ylim([tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale])  
            ax.plot([0,0], [tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
            ax.plot([tdos0[:,0].min()*yscale,tdos0[:,0].max()*yscale], [0,0], color=color[0], linewidth=1.0, alpha=alpha) 
        else:
            plt.ylim([0, tdos0[:,1].max()*yscale])
            ax.plot([0,0], [0, tdos0[:,1].max()*yscale], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
        plt.xlim(elim)
        plt.yticks([])
        
    elif style == 'vertical':
        ax.plot(tdos0[:,1], tdos0[:,0], color=color[0], linewidth=1.1, label='TDOS')
        if pdos0 is not None:
            for orb in range(pdos0.shape[1]): 
                ax.plot(pdos0[:,orb], tdos0[:,0], color=color[orb+1], linewidth=1.0,label=legend[orb])
                if fill == True: ax.fill_between(pdos0[:,orb], tdos0[:,0], color=color[orb+1], alpha=alpha)
            
        if spin == 'updown':
            ax.plot(tdos1[:,1], tdos0[:,0], color=color[0], linewidth=1.1)
            if pdos1 is not None:
                for orb in range(pdos1.shape[1]): 
                    ax.plot(pdos1[:,orb], tdos0[:,0], color=color[orb+1],linewidth=1.0)
                    if fill == True: ax.fill_between(pdos1[:,orb], tdos0[:,0], color=color[orb+1], alpha=alpha)

    
        # Graph adjustments 
        plt.xlabel('DOS', size=fontsize+4)   
        plt.ylabel('Energy (eV)', size=fontsize+4)
        if spin == 'updown':
            plt.xlim([tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale])  
            ax.plot([tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale], [0,0], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
            ax.plot([0,0], [tdos0[:,0].min()*yscale,tdos0[:,0].max()*yscale], color=color[0], linewidth=1.0, alpha=alpha) 
        else:
            plt.xlim([0,tdos0[:,1].max()*yscale])
            ax.plot([0, tdos0[:,1].max()*yscale], [0,0], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
        plt.ylim(elim)
        plt.xticks([])
        
    else:   
        assert 0, "Style must be 'horizontal' or 2 'vertical'"
    # Legend
    lgnd = ax.legend(loc=loc, numpoints=1, fontsize=fontsize)
            
    # Graph adjustments 
    border = 1.08        
    ax.tick_params(labelsize=fontsize, width=border)
    ax.spines['top'].set_linewidth(border)
    ax.spines['right'].set_linewidth(border)
    ax.spines['bottom'].set_linewidth(border)
    ax.spines['left'].set_linewidth(border)
    plt.tight_layout()
    if save == True: 
        fig.savefig(figname+'.'+format, dpi=dpi, format=format)      
    else:
        plt.show() 
        

def plot_kdos(calculator, efermi=None, spin=0, lm=None, plot_band=False, label=None, cmap='Blues', save=False, band_color=['#ffffff','#f2f2f2','#f2f2f2'], figsize=(6,6), figname='BAND', xlim=None, ylim=[-6,6], fontsize=18, dpi=300, format='png'):
    '''Plot k-resolved DOS
       
        Attribute:
            calculator      : the engine to generate band structure (VASP, QE, CRYSTAL, etc.)
            efermi          : a Fermi level or a list of Fermi levels
            spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                              0 or 1 for spin polarized
            color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                              
                              
    '''

    assert isinstance(band_color,list)
    assert len(band_color) == 3
    
    if lm is not None:
        assert isinstance(lm, str) or isinstance(lm, list)
        proj_kdos = True
    else:
        lm = [atom+':s,p,d' for atom in calculator.element]  
        proj_kdos = False
        
    tdos, pdos, kpath, sym_kpoint_coor, formatted_label = calculator._generate_kdos(efermi=efermi, spin=spin, lm=lm, label=label)
    
    if proj_kdos is True:
        kdos = pdos
    else:
        kdos = tdos[:,:,1]
        
    ##----------------------------------------------------------
    ##Plotting:        
    ##----------------------------------------------------------
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    yrange = (-500,500)
               
    # Plot the high symmetric kpoint grid
    if sym_kpoint_coor is not None:
        for kpt in range(sym_kpoint_coor.shape[0]):
            ax.plot([sym_kpoint_coor[kpt]]*2,yrange, color=band_color[1], linewidth=1.0, dashes=[6,3])
                
    if formatted_label is not None and xlim is None:
        nkpts = len(formatted_label)
        assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
        for kpt in range(nkpts):   
            point = formatted_label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(sym_kpoint_coor[kpt]/kpath.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)    

    # Get X, Y, Z
    kpath = np.asarray(kpath).flatten()
    X, Y = np.meshgrid(kpath, tdos[0,:,0])
    Z = kdos.T
    vmin = Z.min()
    vmax = Z.max()

    # Plot kDOS           
    ax.plot([0, kpath.max()], [0,0], color=band_color[2], linewidth=1.0, dashes=[6,3])       # Fermi level
    CS = ax.contourf(X, Y, Z, 100, cmap=cmap, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)) 
    cbar = fig.colorbar(CS, format='%.2f')
    cbar.ax.tick_params(labelsize=fontsize-4)
    
    # Plot band
    if plot_band is True:
        band, kpath, sym_kpoint_coor, formatted_label, conventional = calculator._generate_band(efermi=efermi, spin=spin, label=label)  
        # Plot bands            
        ax.plot([0,kpath.max()],[0,0],color=band_color[2],linewidth=1.0, dashes=[6,3])       # Fermi level
        for ith in range(band.shape[1]):
            ax.plot(kpath.T,band[:,ith],color=band_color[0],linewidth=1.0)    
    
    # Graph adjustments             
    ax.tick_params(labelsize=fontsize)
    if xlim is None:
        plt.xlim([0,kpath.max()])
        plt.xticks([])
        plt.xlabel('k', size=fontsize+4)
    else:
        plt.xlim(xlim)
        plt.xlabel('k ' + r'($\AA^{-1}$)', size=fontsize+4)
    ax.xaxis.set_label_coords(0.5, -0.08) 
    plt.ylabel('Energy (eV)', size=fontsize+4)        
    plt.ylim(ylim)
    plt.tight_layout()
    if save is True: 
        fig.savefig(figname+'.'+format,dpi=dpi,format=format)      
    else:
        plt.show()  
        
