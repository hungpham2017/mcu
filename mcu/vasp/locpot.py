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
from . import utils, vasp_io
import matplotlib as mpl
import matplotlib.pyplot as plt
    
            
class main:
    def __init__(self, locpot='LOCPOT'):
        '''Get LOCPOT file and return a LOCPOT object '''
        self.locpot = vasp_io.LOCPOT(locpot)
        
    def get_2D_average(self, axis='z'):
        return self.locpot.get_2D_average(axis)
        
    def get_vacumm(self, pot=None, axis='z', error=0.01):
        '''Get the electrostatic potential at vacuum
        '''
        if not isinstance(pot,np.ndarray): 
            pot = self.get_2D_average(axis)
        lower_bound = pot[1].max()- 2*error
        idx = pot[1] > lower_bound
        pot_in_window = pot[1,idx]
        e_vacuum = np.average(pot[1,idx])
    
        return e_vacuum
        
    def plot(self, axis='z', error=0.01, color=['r', '#737373'], ylim=None, save=False, figname='elecpot', figsize=(8,6), fontsize=18, dpi=600, format='png'):
        '''Function to plot the inplane average potential to check the convegence
        '''
        
        pot = self.get_2D_average(axis)
        e_vacuum = self.get_vacumm(pot, error=error)
        
        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.plot(pot[0], pot[1], color=color[0],linewidth=1.1,label='TDOS')
        ax.plot([pot[0].min(), pot[0].max()], [e_vacuum,e_vacuum], color=color[1], linewidth=1.1, dashes=[6,3])
        if ylim == None: ylim = (round(pot[1].min()) - 10.0,round(pot[1].max()) + 10.0)
        y_evacuum = (e_vacuum-ylim[0])/(ylim[1]-ylim[0]) 
        label = r'$E_{Vacuum}$ = ' + str(round(e_vacuum,2)) + r' $\pm$ ' + str(round(error,2)) + ' eV' 
        ax.text(0.02, y_evacuum, label  , verticalalignment='bottom', horizontalalignment='left',transform=ax.transAxes, color='black', fontsize=fontsize+4)  
 
        # Graph adjustments
        border = 1.08
        ax.tick_params(labelsize=fontsize, width=border)
        ax.spines['top'].set_linewidth(border)
        ax.spines['right'].set_linewidth(border)
        ax.spines['bottom'].set_linewidth(border)
        ax.spines['left'].set_linewidth(border)
        plt.xlabel(axis + r' ($\AA$)', size=fontsize+4)
        ax.xaxis.set_label_coords(0.5, -0.09) 
        plt.ylabel('Electrostatic potential (V)', size=fontsize+4)        
        plt.ylim(ylim) 
        plt.xlim([0,pot[0].max()])
        plt.tight_layout()
        if save == True: 
            fig.savefig(figname+'.'+format, dpi=dpi, format=format)      
        else:
            plt.show() 
