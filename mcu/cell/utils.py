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
    '''Convert atom from itsZ to its symbol and vice versa''' 
    
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
            for axis in ['x','y','z']: # loop over x, y, z
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
               

def atoms_irred2full(lattice, irred_symbol, irred_frac_coors, rotations, translations, prec=1e-4):
    '''Giving a list of irreducible atoms, generate all the other atoms
    Input:
    =====
        - lattice matrix
        - irreducible atoms symbol
        - irreducible atoms fractional coordinates 
        - Rotation and translation operators in matrix form
        - prec in Angstrom
    '''

    full_atom_Z = []
    full_frac_coors = []
    ntrans = len(translations)
    irred_atom_Z = convert_atomtype(irred_symbol)
    
    # Operate R and T on the irreducible atoms to get a redudant list of atoms
    for i, symbol in enumerate(irred_atom_Z):
        new_atoms = np.einsum('iab,b->ia', rotations, irred_frac_coors[i]) + translations
        full_atom_Z.extend([symbol]*ntrans)
        full_frac_coors.extend(new_atoms)

    # After the symmetry operation, the coordination may be outside the unit cell, i.e., <0 or > 1
    # Here all the atoms are transformed to the unit cell
    natoms = len(full_frac_coors)
    full_frac_coors_inside_cell = np.asarray(full_frac_coors) - np.int64(full_frac_coors)           # Transform >1 to [0,1]
    full_frac_coors = full_frac_coors_inside_cell.flatten()
    negative_coords_idx = np.where(full_frac_coors<0)[0]                                                  # locate negative coordinates
    full_frac_coors[negative_coords_idx] = full_frac_coors[negative_coords_idx] + 1.0                           # Transform <0 to [0,1]
    full_frac_coors = full_frac_coors.reshape(natoms,3)

    # Remove redundant atoms BUT leave the disosdered ones with a Warning for the user
    atoms = np.empty([natoms,4])
    atoms[:,0] = np.asarray(full_atom_Z)
    atoms[:,1:] = full_frac_coors @ lattice             # Use absolute coordinates to check the redudancy
    atoms = np.unique(atoms.round(int(abs(np.log10(prec)))), axis=0)
    
    # Check disodered atoms:
    num_irredundant_atom = atoms.shape[0]
    coords = np.unique(atoms[:,1:].round(int(abs(np.log10(prec)))), axis=0)
    if coords.shape[0] < num_irredundant_atom:
        print("WARNING! Your CIF may contain disordered atoms (Atoms that are too close to each other)")
        return None,  None
    else:    
        irredundant_atom_symbol = convert_atomtype(np.int64(atoms[:,0]))
        irredundant_abs_coors = atoms[:,1:]
        irredundant_frac_coors =  irredundant_abs_coors @ np.linalg.inv(lattice)

        return irredundant_atom_symbol, irredundant_frac_coors
    
    
def convert_lattice(lattice):
    '''Convert a lattice matrix to lattice parameters (a, b, c, alpha, beta, gamma) and vice versa'''
    lattice = np.asarray(lattice)
    if lattice.shape[0] == 3:           # lattice matrix to lattice parameters
        a, b, c = np.linalg.norm(lattice, axis=1)
        cos_alpha = lattice[1].dot(lattice[2])/np.linalg.norm(lattice[1])/np.linalg.norm(lattice[2])
        alpha = np.arccos(cos_alpha)*180/np.pi
        cos_beta = lattice[0].dot(lattice[2])/np.linalg.norm(lattice[0])/np.linalg.norm(lattice[2]) 
        beta = np.arccos(cos_beta)*180/np.pi 
        cos_gamma = lattice[0].dot(lattice[1])/np.linalg.norm(lattice[0])/np.linalg.norm(lattice[1])
        gamma = np.arccos(cos_gamma)*180/np.pi
        
        return np.asarray([a,b,c,alpha,beta,gamma])
        
    elif lattice.shape[0] == 6:         # lattice parameters to lattice matrix
        a,b,c,alpha,beta,gamma = lattice
        vec1 = [a, 0.0, 0.0]    # a aligns along x
        vec2 = [b*np.cos(gamma*np.pi/180), b*np.sin(gamma*np.pi/180), 0.0]    # b is in the xy plane  
        vec3 = [0.0, 0.0, 0.0]
        vec3[0] = np.cos(beta*np.pi/180)*c
        vec3[1] = (np.cos(alpha*np.pi/180)*b*c - vec2[0]*vec3[0])/vec2[1]
        vec3[2] = np.sqrt(c**2 - vec3[0]**2 - vec3[1]**2)  

        return np.asarray([vec1,vec2,vec3])
        

    