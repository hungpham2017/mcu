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

'''
Utilities for vasp module
'''

import os
import numpy as np
from ..utils.misc import check_exist
from ..cell import parameters
from ..cell import utils as cell_utils             
            
    
def str_extract(string, start, end):
    '''Get substring between start and end keyword from a string'''
    
    if not isinstance(start,str) or not isinstance(end,str):
        raise Exception('start and end are two string')
    if  start not in string:
        raise Exception('Cannot find', start, 'in the string:', string)
    if  end not in string:
        raise Exception('Cannot find', end, 'in the string:', string)      
        
    len_start = len(start)
    len_end = len(end)
    
    substring = []    
    stop = False
    new_string = string
    while stop == False:   
        new_string = new_string[(new_string.find(start)+len_start):]
        sub_str = new_string[:new_string.find(end)]
        new_string = new_string[new_string.find(end)+len_end:]
        
        find_start = new_string.find(start)
        find_end = new_string.find(end)        
        if find_start == -1 or find_end == -1: 
            stop = True
        substring.append(sub_str)
       
    if len(substring) == 1 : substring = substring[0] 
    
    return substring
    
def read_WAVEDER(file='WAVEDER'):
    '''Read the WAVEDER file
    
        the matrix of the derivative of the cell periodic part
        of the wavefunctions with respect to i k is:      
        cder = CDER_BETWEEN_STATES(m,n,1:NK,1:ISP,j)= <u_m|  -i d/dk_j | u_n> = - <u_m| r_j |u_n>
    '''
    if not check_exist(file):
        print('Cannot find the %s file. Check the path:' % file)
        
    from scipy.io import FortranFile
    data = FortranFile(file, 'r')
    nb_tot, nbands_cder, nkpts, ispin = data.read_record(dtype= np.int32)
    nodesn_i_dielectric_function = data.read_record(dtype= np.float)
    wplasmon = data.read_record(dtype= np.float).reshape(3,3)
    cder = data.read_record(dtype= np.complex64).reshape(ispin,nkpts,nbands_cder,nb_tot,3)
    
    return cder, nodesn_i_dielectric_function, wplasmon
    
    
def read_WAVEDERF(file='WAVEDERF'):
    '''Read the WAVEDERF file, a formatted WAVEDER file'''
    
    if not check_exist(file):
        print('Cannot find the %s file. Check the path:' % file)
        
    data = open(file, "r").readlines() 
    ispin, nkpts, nbands_cder = np.int32(data[0].split())
    
    # the last index of cder for cdum_x,cdum_y,cdum_z
    cder = np.empty([ispin,nkpts,nbands_cder,nbands_cder,3]) 
    
    line = 1
    for spin in range(ispin):
        for kpt in range(nkpts):
            for band1 in range(nbands_cder):      
                for band2 in range(nbands_cder):      
                    x_real, x_imag, y_real, y_imag, z_real, z_imag = np.float64(data[line].split())[-6:]
                    cdum[spin,kpt,band1,band2] = np.asarray([np.complex(x_real,x_imag), np.complex(y_real,y_imag), np.complex(z_real,z_imag)])     
                    line += 1
                
    return cder
    
def rm_redundant_band(kpts, band):
    '''Remove redundant kpoints from the band calculations'''
    
    redundant_idx = []
    for kpt in range(kpts.shape[0]-1):
        if abs(kpts[kpt] - kpts[kpt+1]).sum() < 1.e-10: redundant_idx.append(kpt+1)
    
    return np.delete(kpts,redundant_idx,axis=0), np.delete(band,redundant_idx,axis=0) 
        
def convert_kpath(kpath):
    '''Provide a kpath string, return a list
    '''

    kpath = kpath.split()
    assert np.mod(len(kpath)-1,4) == 0
    path = kpath[-1].split('-')
    #nkpoint = len(path)
    coor = kpath[:-1]
    ncoor = len(coor)//4
    highk_data = [[coor[i*4], np.float64(coor[i*4+1:i*4+4])] for i in range(ncoor)]
    k_list = []    
    for k in path:
        find = False
        for highk in highk_data:
            if k == highk[0]: 
                k_list.append([k, highk[1]])
                find = True
        if find == False: raise Exception("Cannot find", k, "in the high symmetric kpoint provided")
                
    return k_list
   
def get_1Dkpath(kpath, npoint=20):
    '''Provide high symmetric kpoint coordinates and labels and a path, this function returns a KPOINT file for band structure computation'''
    
    k_list = convert_kpath(kpath)
    
    temp = np.arange(npoint+1).reshape(-1,1)
    kpts = []
    npath = len(k_list) - 1
    for path in range(1, len(k_list)):
        line = k_list[path][1] - k_list[path-1][1]
        kpts.append(k_list[path-1][1] + temp*line/npath)
    kpts = np.asarray(kpts).reshape(npath*(npoint+1),-1)
    kpts = rm_redundant_band(kpts,kpts)[0]
    
    with open('KPOINTS', 'w') as f:
        f.write('Generated mesh by mcu\n')	
        f.write('   %4d + # of k-points from IBZKPT\n' % (kpts.shape[0]))	    
        f.write('Reciprocal lattice\n')
        f.write('   #This is where you add the k-mesh copied from IBZKPT file, this k-mesh is used to run SCF\n')
        for k in range(kpts.shape[0]):
            f.write('%10.7f  %10.7f  %10.7f %2d\n' % (kpts[k,0], kpts[k,1], kpts[k,2], 0))
    
    return kpts
    
def cell_to_spgcell(cell, atoms):
    '''Providing the cell attribute from vasprun.xml, return the cell for spglib'''
    
    lattice = cell[0]
    positions = cell[2]
    numbers = cell_utils.convert_atomtype(atoms)
    return (lattice, positions, numbers)
