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
Utilities for cell
'''

import os, datetime
import numpy as np
from mcu.cell import parameters            


def to_symbol(atoms_number):
    '''Convert atom from Z to symbol, both are lists''' 

    atom = []
    for Z in atoms_number:
        stop = False
        for element in parameters.ELEMENTS.items():
            if Z == element[1][0]:
                atom.append(element[0])
                break
    return atom
    
def convert_frac(frac, err=1.e-5):
    '''Convert a fraction (str) to float and vice versa'''
    
    if isinstance(frac,str):
        frac = frac.split("/")
        if len(frac) == 1: 
            frac = float(frac)
        else:
            frac = float(frac[0]) / float(frac[1])
    else:
        
        if abs(frac - 1/2) < err: 
            frac = '1/2'
        elif abs(frac - 1/3) < err: 
            frac = '1/3'
        elif abs(frac - 2/3) < err: 
            frac = '2/3'
        elif abs(frac - 1/6) < err: 
            frac = '1/6' 
        elif abs(frac - 5/6) < err: 
            frac = '5/6'
        elif (abs(frac) < err) or abs(frac - 1) < err: 
            frac = '0'
        
    return frac
        
def symop_xyz2mat(sym_operators):
    '''Convert string operator to the matrix form'''
    
    rotations = []
    translations = []
    
    
    for opt in sym_operators: # loop over each sym opt
        opt = opt.split(",")
        rotation = []
        translation = []
        for i in range(3):  # loop over rows of rotation mat
            temp = opt[i]
            vec = []
            axes = ['x','y','z']
            for axis in axes: # loop over axis of row vector
                if axis in temp:
                    idx = temp.index(axis)
                    if idx == 0:
                        sign = +1
                        temp = temp.replace(axis,'')
                    else:
                        sign = temp[idx-1]
                        temp = temp.replace(sign+axis,'')
                        if sign == '+': sign = +1
                        if sign == '-': sign = -1
                else:
                    sign = 0
                vec.append(sign)
            if temp == '':
                temp = 0
            else:
                temp = convert_frac(temp)
            translation.append(temp)
            rotation.append(vec)
        rotations.append(rotation)
        translations.append(translation) 
        
    return rotations, translations
           
def symop_mat2xyz(rotations, translations):
    '''Convert mat operator to the string form'''
    
    num_opt = len(translations)
    syms = []    
    for opt in range(num_opt):   # loop over each sym opt
        rotation = rotations[opt]
        translation = translations[opt]        
        sym = ''
        for i in range(3):      # loop over rows of rotation mat
            # Adding the rotation opt
            axes = ['x','y','z']
            for j in range(3): # loop over axis of row vector
                temp1 = rotation[i][j]
                if temp1 < 0:
                    sym += '-' + axes[j]
                elif temp1 > 0:
                    sym += '+' + axes[j] 
            # Adding the translation opt   
            tran = convert_frac(abs(translation[i]))
            if translation[i] < 0: sign = '-'
            if translation[i] > 0: sign = '+'
            if abs(translation[i]) < 1.e-8: sign = ''
            if tran != '0':
                sym = sym + sign + tran
            sym += ","
        syms.append(sym[:-1])
        
    return syms
                    
def rm_paren(string):
    return string.replace("(","").replace(")","")
    
def convert_lattice(lattice):
    '''Convert a lattice matrix to cif file format cell'''    
    lattice = np.asarray(lattice)
    a, b, c = np.linalg.norm(lattice, axis=1)
    cos_alpha = lattice[1].dot(lattice[2])/np.linalg.norm(lattice[1])/np.linalg.norm(lattice[2])
    alpha = np.arccos(cos_alpha)*180/np.pi
    cos_beta = lattice[0].dot(lattice[2])/np.linalg.norm(lattice[0])/np.linalg.norm(lattice[2]) 
    beta = np.arccos(cos_beta)*180/np.pi 
    cos_gamma = lattice[0].dot(lattice[1])/np.linalg.norm(lattice[0])/np.linalg.norm(lattice[1])
    gamma = np.arccos(cos_gamma)*180/np.pi
    
    return (a,b,c,alpha,beta,gamma)
    
def frac2cart(cell):
    '''Convert fractional coordinates to cartesian
        Attributes:
            cell        : a celll object in spglib format
        Return:
            postions    : in cartesian coordinates
    '''
    
    lattice = cell[0]    
    frac_coordinates = cell[1]
    lattice_length = np.sqrt((np.sum(lattice**2, axis = 1)))
    cart_coordinates = frac_coordinates*lattice_length
    return cart_coordinates
    