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
from matplotlib.ticker import AutoMinorLocator

def analyze(elastic_tensors):
    ''' Ref: 10.1021/jacs.8b13075 '''
    print("Unit: N/m") 
    print("  C11     C22    C12     C66     Gamma    Y_x     Y_y     v_x     v_y")  
    for tensor in elastic_tensors:
        C11, C22, C12, C66 = tensor
        C21 = C12
        # Layer modulus
        gamma = 1./4.*(C11 + C22 + 2*C12)

        # Young modulus
        Y_x = (C11*C22 - C12*C21) / C22
        Y_y = (C11*C22 - C12*C21) / C11

        # Poisson's ratio
        v_x = C21/C22
        v_y = C12/C11
        
        print("{0:3.2f}   {1:3.2f}  {2:3.2f}  {3:3.2f}  {4:3.2f}  {5:3.2f}  {6:3.2f}  {7:2.4f}  {8:2.4f}".format(C11, C22, C12, C66, gamma, Y_x, Y_y, v_x, v_y))
    
# In-plane Young's modules Y(theta), theta: 0 -> 360 degree, θ is the angle relative to the x direction
def young_theta(theta, elastic_tensor):
    C11, C22, C12, C66 = elastic_tensor
    numerator = C11*C22 - C12**2
    denominator = C11*(np.sin(theta))**4 + C22*(np.cos(theta))**4 \
        + (numerator/C66 - 2*C12) * (np.cos(theta))**2 * (np.sin(theta))**2
    return numerator/denominator
    

# Poisson's ratio v(theta)
def poisson_theta(theta, elastic_tensor):
    C11, C22, C12, C66 = elastic_tensor
    a = C11*C22 - C12**2
    numerator = (C11 + C22 - a/C66) * (np.cos(theta))**2 * (np.sin(theta))**2 \
                - C12 * ((np.sin(theta))**4 + (np.cos(theta))**4)
    
    denominator = C11*(np.sin(theta))**4 + C22*(np.cos(theta))**4 \
        + (a/C66 - 2*C12) * (np.cos(theta))**2 * (np.sin(theta))**2
    return -numerator/denominator
    

# Plotting polar diagram for in-plane Young's modules
def plot_polar(elastic_tensors, young=True, legend=None, save=False, thetagrids=3, color=None, 
                gridcolor='#a6a6a6', grid_style=':', figsize=5.0, figname='polardiagram', 
                fontsize=12, dpi=300, format='png'):
                              
    '''       
    Attributes:
        elastic_tensors : a list of elastic tensor (4 components vectors)
        thetagrids      : 1, 2, 3, ...
        poisson         : True if Poisson's ratio is plotted instead, othewise plotting Young's modules
    
    '''
    if color is None:
        color = ['r','g','b','y','m','c']
        
    # Estimate "best" parameters:
    letter_size = fontsize/72
    shift_left = 6.0*letter_size
    fig_h = figsize
    fig_w = 0.8*fig_h + 2.5*shift_left
    
    if legend is not None:
        assert len(legend) == len(elastic_tensors), "The number of legends must be equal to the number of elastic tensor"
        longest_string = len(max(legend, key=len))
        fig_w += longest_string*letter_size

    polar_shift = (1.5*shift_left + 0.8*fig_h/2.0 - 0.8*fig_w/2.0) / fig_w
    theta = np.arange(0,360,1.0) * np.pi / 180.0
    fig = plt.figure(figsize=(fig_w,fig_h))     #5*letter_size/fig_w
    ax = fig.add_subplot(111, projection='polar', position=(polar_shift,.1,0.8,0.8))
    r_max = 0.0
    func = young_theta
    if young == False: func = poisson_theta
    for i, tensor in enumerate(elastic_tensors):
        r = func(theta, tensor)
        if legend is not None:
            label = legend[i]
        else:
            label = None
        ax.plot(theta, r, color=color[i], label=label)
        if r.max() > r_max: r_max = r.max()

    if young:
        r_max = 50.0 * np.ceil(r_max/50.0)
        ax.set_rmax(r_max)
        ax.set_rticks(np.arange(1,round(r_max/50.0)+1,1)*50.0)  # Less radial ticks
        
        border =1.5
        ax.spines['polar'].set_linewidth(border)
        
        # Add Cartesian axes
        ax2 = fig.add_axes((shift_left/fig_w,0.1,0.0,0.8))
        ax2.xaxis.set_visible(False) # hide x axis
        nstick = int(round(r_max/100.0))
        ratio = 100.0*nstick/r_max
        yticks_location = np.linspace(0, 1, 2*nstick + 1)*ratio + (1 - ratio)*0.5
        ax2.set_yticks(yticks_location) # set new tick positions 
        ax2.set_yticklabels(list(map(lambda x: str(x), np.arange(- nstick - 1, nstick + 1, 1)[1:]*100)), fontsize=fontsize)
        ax2.yaxis.set_minor_locator(AutoMinorLocator(2)) # set minor tick for every second tick
        plt.ylabel("Y" + r"($\theta$)" + " (N/m)", size=1.2*fontsize)  
    else:
        r_max = 100*r_max
        r_max = 5.0 * np.ceil(r_max/5.0)
        ax.set_rmax(r_max/100.0)
        ax.set_rticks(np.arange(1,round(r_max/5.0)+1,1)*5.0/100.0)  # Less radial ticks
        
        # Add Cartesian axes
        ax2 = fig.add_axes((shift_left/fig_w,0.1,0.0,0.8))
        ax2.xaxis.set_visible(False) # hide x axis
        nstick = int(round(r_max/10.0))
        ratio = 10.0*nstick/r_max
        yticks_location = np.linspace(0, 1, 2*nstick + 1)*ratio + (1 - ratio)*0.5
        ax2.set_yticks(yticks_location) # set new tick positions 
        ax2.set_yticklabels(list(map(lambda x: str(x), np.arange(- nstick - 1, nstick + 1, 1)[1:] / 10.0)), fontsize=fontsize)
        ax2.yaxis.set_minor_locator(AutoMinorLocator(2)) # set minor tick for every second tick
        plt.ylabel(r"$\nu(\theta)$" + " (N/m)", size=1.2*fontsize)  
    
    thetagrids = np.int64(np.arange(360.0/(90.0/thetagrids))*(90.0/thetagrids))
    ax.set_thetagrids(angles=thetagrids, labels=list(map(lambda x: str(x) + u'\N{DEGREE SIGN}', thetagrids)), fontsize=fontsize)
    ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line 
    ax.grid(True, color=gridcolor, linestyle=grid_style)
    ax.set_yticklabels([])
    
    if legend is not None:
        shift_x = (fig_w - 2.0*shift_left) / 0.8 / fig_h
        ax.legend(loc="upper right", bbox_to_anchor=(shift_x, 1.0), numpoints=1, fontsize=fontsize)
    
    if save == True: 
        fig.savefig(figname +'.'+ format,dpi=dpi,format=format)      
    else:
        plt.show()
