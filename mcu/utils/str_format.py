#!/usr/bin/env python -u 
'''
pDMET: Density Matrix Embedding theory for Periodic Systems
Copyright (C) 2018 Hung Q. Pham. All Rights Reserved.
A few functions in pDMET are modifed from QC-DMET Copyright (C) 2015 Sebastian Wouters

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
from itertools import permutations

'''Format orbitals '''

def generate_lm_dict(type='lm'):
    if type == 'lm':
        lm_basis_dict = {'s':['s'], 'p':['px','py','pz'], 'd':['dxy', 'dyz','dz2','dxz','dx2-y2'],  
                         'f':['fz3', 'fxz2', 'fyz2', 'fz(x2-y2)', 'fxyz', 'fx(x2-3y2)', 'fy(3x2-y2)']}
        lm_basis_list = sum(list(lm_basis_dict.values()),[])
        lm_shortcut = []
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 1)])
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 2)])
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 3)])
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 4)])
        lm_shortcut = sum(lm_shortcut,[])
    elif type == 'l':
        lm_basis_dict = {'s':['s'], 'p':['p'], 'd':['d'],  'f':['f']}
        lm_basis_list = sum(list(lm_basis_dict.values()),[])
        lm_shortcut = []
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 1)])
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 2)])
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 3)])
        lm_shortcut.append([''.join(i) for i in permutations(['s', 'p', 'd', 'f'], 4)])
        lm_shortcut = sum(lm_shortcut,[])
    elif type == 'ud':
        lm_basis_dict = {'u':['u'], 'd':['d']}
        lm_basis_list = sum(list(lm_basis_dict.values()),[])
        lm_shortcut = ['ud', 'du']
        
    elif type == 'spin':
        lm_basis_dict = {'s':['s']}
        lm_basis_list = sum(list(lm_basis_dict.values()),[])
        lm_shortcut = ['s']
        
    return lm_basis_dict, lm_basis_list, lm_shortcut

def basic_lm(lm, type='lm'):
    lm_basis_dict, lm_basis_list, lm_shortcut = generate_lm_dict(type)

    formatted_lm = []
    if lm in lm_basis_list:
        formatted_lm.append([lm])
    elif lm in lm_shortcut:
        for orb in lm:
            formatted_lm.append(lm_basis_dict[orb])
    else:
        assert 0, "Cannot recognize the following orbital: " + lm
        
    formatted_lm = sum(formatted_lm,[])

    return formatted_lm
        
def str_lm(lm, type='lm'):
    '''Giving a lm string, return a formatted list which consists of sub lists.
       The quantity (band or DOS) will be summed over all the orbitals in each sub-list 
    '''  
    assert isinstance(lm,str), "lm must be a string"
    lm_basis_dict, lm_basis_list, lm_shortcut = generate_lm_dict(type)
    
    lm = lm.strip()
    
    if lm in lm_basis_list:         # e.g.: lm = 'dx2-y2' or lm = 'p'
        formatted_lm = [basic_lm(lm, type)]
        formatted_atom = [None]
        
    elif lm in lm_shortcut:         
        formatted_lm = [basic_lm(lm, type)]
        formatted_atom = [None]
        
    elif ';' in lm:                 # e.g.: lm = ' 'C:s,p ; Si:sp' and space doesn't matter
        lm_group  = lm.split(';') 
        formatted_lm = []
        formatted_atom = []
        for i, each_lm in enumerate(lm_group):     # e.g.: 'C:s,p'
            each_lm = each_lm.strip()
            if each_lm in lm_basis_list: 
                atom = None
                formatted_lm.append(basic_lm(each_lm, type))
            elif each_lm in lm_shortcut:
                atom = None
                formatted_lm.append(basic_lm(each_lm, type))
            elif each_lm.count(':') == 1: 
                atom, lm_ = each_lm.split(':')
                atom = atom.strip()
                lm_ = lm_.strip()
                formatted_lm_ = []
                if ',' in lm_:
                    temp  = lm_.split(',') 
                    for i, orbs in enumerate(temp):   
                        orbs = orbs.strip()
                        formatted_lm_.append(basic_lm(orbs, type))
                else:
                    formatted_lm_.append(basic_lm(lm_, type))
                formatted_lm_ = sum(formatted_lm_,[])
                formatted_lm.append(list(dict.fromkeys(formatted_lm_)))
            elif (',') in each_lm:
                temp  = each_lm.split(',') 
                formatted_lm_ = []
                for i, orbs in enumerate(temp):
                    orbs = orbs.strip()
                    formatted_lm_.append(basic_lm(orbs, type))
                formatted_lm_ = sum(formatted_lm_,[])
                formatted_lm.append(list(dict.fromkeys(formatted_lm_)))        
                atom = None
            else:
                atom = each_lm
                formatted_lm.append([None])
            formatted_atom.append(atom)
            
    elif lm.count(':') == 1:
        atom, lm_ = lm.split(':')
        atom = atom.strip()
        lm_ = lm_.strip()
        formatted_lm = []
        formatted_atom = [atom]
        if ',' in lm_:
            temp  = lm_.split(',') 
            for i, orbs in enumerate(temp):  
                orbs = orbs.strip()
                formatted_lm.append(basic_lm(orbs, type))
        else:
            formatted_lm.append(basic_lm(lm_, type))     
        formatted_lm = sum(formatted_lm,[])
        formatted_lm = [list(dict.fromkeys(formatted_lm))]    # remove redundant orbitals
        
    elif (',') in lm:
        temp  = lm.split(',') 
        formatted_atom = [None]
        formatted_lm_ = []
        for i, orbs in enumerate(temp): 
            orbs = orbs.strip()
            formatted_lm_.append(basic_lm(orbs, type))
        formatted_lm_ = sum(formatted_lm_,[])
        formatted_lm = [formatted_lm_]

    else:
        formatted_atom = [lm]
        formatted_lm = [[None]]
        
    return formatted_atom, formatted_lm 

def list_lm(lms, type='lm'):
    '''Giving a list of lm strings, return a formatted list which consists of sub lists.
       The quantity (band or DOS) will be summed over all the orbitals in each sub-list  
    '''  
    assert isinstance(lms,list), "lm must be a list of lm strings"
    formatted_atom = []
    formatted_lm = [] 
    for i, lm in enumerate(lms):
        assert isinstance(lm,str), lm + " is not a lm string"
        atom_, formatted_lm_ = str_lm(lm, type)
        formatted_atom.append(atom_)
        formatted_lm.append(formatted_lm_)
        
    return formatted_atom, formatted_lm

def general_lm(lm, type='lm'):
    '''Giving a string or a list of lm strings, return a formatted list which consists of sub lists.
       The quantity (band or DOS) will be summed over all the orbitals in each sub-list  
    '''  
    if isinstance(lm, str):
        formatted_atom = []
        formatted_lm = []
        formatted_atom_, formatted_lm_ = str_lm(lm, type)
        for i, each_lm in enumerate(formatted_lm_):
            formatted_lm.append([each_lm])
            formatted_atom.append([formatted_atom_[i]])
        return formatted_atom, formatted_lm
        
    elif isinstance(lm, list): 
        return list_lm(lm, type)
        
    else:
        assert 0, "lm must be a string or a list of strings" 
    
def format_atom(atom):    
    '''Giving a atom string, return the element and the id'''
    elememnt = "".join([letter for letter in atom if not letter.isdigit()])
    id = "".join([letter for letter in atom if letter.isdigit()])    
    if id is '':
        id = None
    else:
        id = int(id)
    return elememnt, id
        
def format_klabel(klabel):
    '''Giving a kpath label string, return the label list
        Each k is separated by ";" or "-"
        '''
    assert isinstance(klabel, str), "Labels for k-point must be a string. Check it!"
    
    klabel = klabel.strip()
    temp = klabel.replace(";"," ").replace("-"," ").split()
    is_digit = sum([item.replace(".","").isdigit() for item in temp])
    if is_digit == 0:
        labels = temp
        coords = None
    else:
        assert not np.mod(len(temp),4) , 'Wrong format for the k-point labels. Check it!'
        nkpt = len(temp) // 4
        labels = [temp[i] for i in np.arange(nkpt)*4]
        coords = np.float64([[temp[i + 1], temp[i + 2], temp[i + 3]] for i in np.arange(nkpt)*4]) 
        assert -1 <= coords.all() <= 1, 'k-point coordinates needs to be in the fractional unit.'
    return labels, coords
    
def format_legend(legend):
    '''Giving a legend string, return the label list'''
    if isinstance(legend, list):
        for item in legend:
            assert isinstance(item, str), "Legend needs to be a list of string or a string"
        return legend
    elif isinstance(legend, str):    
        legend = legend.strip()
        temp = legend.split(";")
        out = [label.strip() for label in temp]
        return out
    else:
        assert 0, "Legend needs to be a list of string or a string"