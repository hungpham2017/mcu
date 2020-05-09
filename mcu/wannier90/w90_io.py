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
import re, textwrap
from ..utils.misc import check_exist
from ..cell import utils as cell_utils  



def get_info_from_block(data, key):
    '''Provide data and a key, return info from the block indicated by key'''
    block_MATCH = re.compile(r'''
    [\w\W]*
      (?i)begin [ ]* ''' + key.strip() + '''
    (?P<content>
      [\s\S]*?(?=\n.*?[ ] $|(?i)end)  # match everything until next blank line or EOL
    )
    ''', re.VERBOSE)
    match = block_MATCH.match(data)
    if match is not None:
        return match['content']
    else:
        return match

def get_info_after_key(data, key):
    '''Provide data and a key, return info from the block indicated by key'''
    key_MATCH = re.compile(r'''
    [\w\W]* ''' + key.strip() + ''' [ ]* [:=]+''' + '''(?P<value>[\s\S]*?(?=\n.*?[ ] $|[\n]))''', re.VERBOSE)
    match = key_MATCH.match(data)
    if match is not None:
        return match['value']
    else:
        return match

def read_win(filename):
    '''Read seedname.win'''

    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, 'r') as data_file:
        data = data_file.read()
        
        unit_cell = get_info_from_block(data, 'Unit_Cell_Cart')
        unit_cell = np.float64(unit_cell.split()).reshape(3,3)
        
        abs_coords = None  
        atoms_cart = get_info_from_block(data, 'atoms_cart')
        if atoms_cart is not None:
            atoms_cart = atoms_cart.split()
            natom = len(atoms_cart) // 4
            atom = [atoms_cart[4*i] for i in range(natom)]
            abs_coords = np.float64([atoms_cart[4*i + 1 : 4*i + 4] for i in range(natom)])
            
        frac_coords = None    
        atoms_frac = get_info_from_block(data, 'atoms_frac')
        if atoms_frac is not None:
            atoms_frac = atoms_frac.split()
            natom = len(atoms_car) // 4
            atom = [atoms_frac[4*i] for i in range(natom)]
            frac_coords = np.float64([atoms_frac[4*i + 1 : 4*i + 4] for i in range(natom)])            
        
        mp_grid = np.int64(get_info_after_key(data, 'mp_grid').split()).tolist()
        
        kpoint_path = get_info_from_block(data, 'kpoint_path')
        kpath = None
        if kpoint_path is not None:
            kpoint_path_MATCH = re.compile(r'''
            [ ]* 
            (?P<k1>\S+) [ ]* (?P<k1x>\S+) [ ]* (?P<k1y>\S+) [ ]* (?P<k1z>\S+) [ ]* (?P<k2>\S+) [ ]* (?P<k2x>\S+) [ ]* (?P<k2y>\S+) [ ]* (?P<k2z>\S+)
            ''', re.VERBOSE)
            kpath = []
            for kpoint in kpoint_path_MATCH.finditer(kpoint_path):
                content = kpoint.groupdict()
                k1 = [content['k1'], np.float64([content['k1x'], content['k1y'], content['k1z']])]
                k2 = [content['k2'], np.float64([content['k2x'], content['k2y'], content['k2z']])]
                kpath.append([k1, k2])
            
        kpts = None 
        kpts_data = get_info_from_block(data, 'kpoints')
        kpts = np.float64(kpts_data.split()).reshape(-1, 3)

        out = {}
        out['unit_cell'] = unit_cell
        out['atom'] = atom
        out['abs_coords'] = abs_coords
        out['frac_coords'] = frac_coords
        out['mp_grid'] = mp_grid
        out['kpath'] = kpath
        out['kpts'] = kpts
        
        return out  


def read_band(filename):
    '''Read seedname_band.dat'''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, 'r') as data_file:
        data = data_file.read()
        temp = data.split('\n  \n')[:-1]
        bands = []
        for i, band in enumerate(temp):
            formatted_band = np.float64(band.split()).reshape(-1, 2)    
            if i == 0: proj_kpath = formatted_band[:,0]
            bands.append(formatted_band[:,1])

        return proj_kpath, np.asarray(bands).T

def read_kpt(filename):
    '''Read seedname_kpt.dat'''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, 'r') as data_file:
        data = data_file.read()
        temp = data.split()
        nkpts = int(temp[0])
        kpath_frac = np.float64(temp[1:]).reshape(-1, 4)
        return kpath_frac
        