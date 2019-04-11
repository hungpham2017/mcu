#!/usr/bin/python
# -*- coding: utf-8 -*-

'''
Utilities
'''

import os, datetime
import numpy as np
            
            
def check_exist(file):
    '''Check if a file exists in the running directory '''        
    exist = os.path.exists(file)
    return exist
    
def str_extract(string, start, end):
    '''Get sbustring between start and end keyword from a string'''
    
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
    
def read_WAVEDER(waveder = 'WAVEDER'):
    '''Read the WAVEDER file
    
        the matrix of the derivative of the cell periodic part
        of the wavefunctions with respect to i k is:      
        cder = CDER_BETWEEN_STATES(m,n,1:NK,1:ISP,j)= <u_m|  -i d/dk_j | u_n> = - <u_m| r_j |u_n>
    '''
    
    from scipy.io import FortranFile
    file = FortranFile(waveder, 'r')
    nb_tot, nbands_cder, nkpts, ispin = file.read_record(dtype= np.int32)
    nodesn_i_dielectric_function = file.read_record(dtype= np.float)
    wplasmon = file.read_record(dtype= np.float).reshape(3,3)
    cder = file.read_record(dtype= np.complex64).reshape(ispin,nkpts,nbands_cder,nbands_cder,3)
    
    return cder, nodesn_i_dielectric_function, wplasmon
    
    
def read_WAVEDERF(wavederf = 'WAVEDERF'):
    '''Read the WAVEDERF file, a formatted WAVEDER file'''
    
    file = open(wavederf, "r").readlines() 
    ispin, nkpts, nbands_cder = np.int32(file[0].split())
    
    # the last index of cder for cdum_x,cdum_y,cdum_z
    cder = np.empty([ispin,nkpts,nbands_cder,nbands_cder,3]) 
    
    line = 1
    for spin in range(ispin):
        for kpt in range(nkpts):
            for band1 in range(nbands_cder):      
                for band2 in range(nbands_cder):      
                    x_real, x_imag, y_real, y_imag, z_real, z_imag = np.float64(file[line].split())[-6:]
                    cdum[spin,kpt,band1,band2] = np.asarray([np.complex(x_real,x_imag), np.complex(y_real,y_imag), np.complex(z_real,z_imag)])     
                    line += 1
                
    return cder
    
def rm_redundant_band(kpts, band):
    '''Remove redundant kpoints from the band calculations'''
    
    redundant_idx = []
    for kpt in range(kpts.shape[0]-1):
        if abs(kpts[kpt] - kpts[kpt+1]).sum() < 1.e-10: redundant_idx.append(kpt+1)
    
    return np.delete(kpts,redundant_idx,axis=0), np.delete(band,redundant_idx,axis=0) 
        
    
    
    
        
    
    


