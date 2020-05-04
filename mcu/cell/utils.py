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
from . import parameters            


def convert_atomtype(Z_or_symbol):
    '''Convert atom from Z to symbol and vice versa''' 
    
    # Z to symbol
    if not isinstance(Z_or_symbol[0], str):
        symbol = []
        for Z in Z_or_symbol:
            stop = False
            for element in parameters.ELEMENTS.items():
                if Z == element[1][0]:
                    symbol.append(element[0])
                    break
        return symbol
        
    # symbol to Z
    else:
        Z = [parameters.ELEMENTS[atom][0] for atom in Z_or_symbol]
        return Z
    
    
def convert_frac(frac, err=1.e-5):
    '''Convert a fraction (str) to float and vice versa'''
    
    if isinstance(frac,str):
        frac = frac.split("/")
        if len(frac) == 1: 
            frac = float(frac[0])
        else:
            frac = float(frac[0]) / float(frac[1])
    else:
        if abs(frac - 1/2) < err: frac = '1/2'
        elif abs(frac - 1/3) < err: frac = '1/3'
        elif abs(frac - 2/3) < err: frac = '2/3'
        elif abs(frac - 1/4) < err: frac = '1/4'
        elif abs(frac - 3/4) < err: frac = '3/4'
        elif abs(frac - 5/4) < err: frac = '5/4'
        elif abs(frac - 1/6) < err: frac = '1/6' 
        elif abs(frac - 5/6) < err: frac = '5/6'
        elif (abs(frac) < err) or abs(frac - 1) < err: frac = '0'
        else:
            frac = False
            
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
            temp = opt[i].strip()
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
            if tran is False:
                print("WARNING! Cannot converge the translations of " + str(translation[i]) + " to the fractional value used by CIF. The resulting symmetry operator in CIF is wrong")
            elif tran != '0':
                sym = sym + sign + tran
            sym += ","
        syms.append(sym[:-1])
        
    return syms
               

def genetate_atoms(irred_symbol, irred_frac, rotations, translations, prec=8):
    '''Operate the R and T operators on irreducible atoms then remove redundant ones'''
    
    full_Z = []
    full_frac = []
    nsymopt = len(translations)
    
    irred_Z = convert_atomtype(irred_symbol)
    # Operate R and T on the irreducible atoms
    for i, Z in enumerate(irred_Z):
        new_atoms = np.einsum('iab,b->ia', rotations, irred_frac[i]) + translations
        full_Z.extend([Z]*nsymopt)
        full_frac.extend(new_atoms)

    natoms = len(full_frac)
    full_frac = np.asarray(full_frac) - np.int64(full_frac)
    full_frac = full_frac.flatten()
    full_frac[np.where(full_frac<0)[0]] = full_frac[np.where(full_frac<0)[0]] + 1.0
    full_frac = full_frac.reshape(natoms,3)
    atoms = np.empty([natoms,4])
    atoms[:,0] = np.asarray(full_Z)
    atoms[:,1:] = full_frac
    atoms = np.unique(atoms.round(prec), axis=0)
    full_Z = np.ndarray.tolist(np.int64(atoms[:,0]))
    full_frac = atoms[:,1:]

    return full_Z, full_frac
    
def rm_paren(string):
    return string.replace("(","").replace(")","")
    
    
def convert_lattice(lattice):
    '''Convert a lattice matrix to cif file format cell and vice versa'''
    
    lattice = np.asarray(lattice)
    if lattice.shape[0] == 3:
        a, b, c = np.linalg.norm(lattice, axis=1)
        cos_alpha = lattice[1].dot(lattice[2])/np.linalg.norm(lattice[1])/np.linalg.norm(lattice[2])
        alpha = np.arccos(cos_alpha)*180/np.pi
        cos_beta = lattice[0].dot(lattice[2])/np.linalg.norm(lattice[0])/np.linalg.norm(lattice[2]) 
        beta = np.arccos(cos_beta)*180/np.pi 
        cos_gamma = lattice[0].dot(lattice[1])/np.linalg.norm(lattice[0])/np.linalg.norm(lattice[1])
        gamma = np.arccos(cos_gamma)*180/np.pi
        
        return np.asarray([a,b,c,alpha,beta,gamma])
        
    elif lattice.shape[0] == 6:
        a,b,c,alpha,beta,gamma = lattice
        vec1 = [a, 0.0, 0.0]    # a aligns along x
        vec2 = [b*np.cos(gamma*np.pi/180), b*np.sin(gamma*np.pi/180), 0.0]    # b is in the xy plane  
        vec3 = [0.0, 0.0, 0.0]
        vec3[0] = np.cos(beta*np.pi/180)*c
        vec3[1] = (np.cos(alpha*np.pi/180)*b*c - vec2[0]*vec3[0])/vec2[1]
        vec3[2] = np.sqrt(c**2 - vec3[0]**2 - vec3[1]**2)  

        return np.asarray([vec1,vec2,vec3])
        

    