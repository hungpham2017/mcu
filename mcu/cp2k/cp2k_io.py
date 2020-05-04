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
import re, argparse
from ..utils.misc import check_exist
from ..cell import utils as cell_utils  

#--------------------------------------------------------------------------------------------------------- 
# USEFUL FUNCTIONS TO MANIPULATE CP2K files
# the read_bs function was modified from cp2k_bs2csv.py
# ref: https://www.cp2k.org/exercises:2016_uzh_cmest:band_structure_calculation
#---------------------------------------------------------------------------------------------------------
def copy_block(blocks, keyword="CELL_TOP"):
    '''Copy a block indetified by the keyword'''


    start_key = keyword
    block = []
    copy_start = False
    copy_stop = False
    for line in blocks:
        if start_key == line.strip()[:len(start_key)]: 
            block.append(line)
            copy_start = True
        elif copy_start == True:
            copy_stop = True
        if copy_stop: break
              
    return block
    
def get_value(blocks, keyword="Fermi energy"):
    '''Get value after a keyword'''


    start_key = keyword
    for line in blocks:
        if start_key == line.strip()[:len(start_key)]: 
            value = np.float64(line.split()[-1])
            break     
    return value 
    
def get_atoms(blocks):
    "Get the atom type and coordinates"
    
    nline = len(blocks)
    start_key = "Atom  Kind  Element"
    atom_type = []
    atom_position = []
    copy_stop = False
    for i, line in enumerate(blocks):
        if start_key == line.strip()[:len(start_key)]: 
            for j in range(i + 2, nline):
                temp = blocks[j].split()
                if temp == []:
                    break
                else:
                    atom_type.append(temp[2])
                    atom_position.append(temp[4:7]) 
            copy_stop = True                    
        elif copy_stop == True:
            break
            
    return atom_type, np.float64(atom_position)
                      
#Used in the read_bs function
SET_MATCH = re.compile(r'''
[ ]*
  SET: [ ]* (?P<setnr>\d+) [ ]*
  TOTAL [ ] POINTS: [ ]* (?P<totalpoints>\d+) [ ]*
  \n
(?P<content>
  [\s\S]*?(?=\n.*?[ ] SET|$)  # match everything until next 'SET' or EOL
)
''', re.VERBOSE)

SPOINTS_MATCH = re.compile(r'''
[ ]*
  POINT [ ]+ (?P<pointnr>\d+) [ ]+ (?P<a>\S+) [ ]+ (?P<b>\S+) [ ]+ (?P<c>\S+)
''', re.VERBOSE)

POINTS_MATCH = re.compile(r'''
[ ]*
  Nr\. [ ]+ (?P<nr>\d+) [ ]+
  Spin [ ]+ (?P<spin>\d+) [ ]+
  K-Point [ ]+ (?P<a>\S+) [ ]+ (?P<b>\S+) [ ]+ (?P<c>\S+) [ ]*
  \n
[ ]* (?P<npoints>\d+) [ ]* \n
(?P<values>
  [\s\S]*?(?=\n.*?[ ] Nr|$)  # match everything until next 'Nr.' or EOL
)
''', re.VERBOSE)

        
def read_ouput(filename):  
    '''Read the cp2k output file'''
    
    assert check_exist(filename), "Cannot find " + filename
    with open(filename, "r") as file:
        outfile = file.read().splitlines()
        
    # Cell information
    cell_block = copy_block(outfile, "CELL_TOP")
    a_vec =cell_block[1].split()[4:7]
    b_vec = cell_block[2].split()[4:7]
    c_vec = cell_block[3].split()[4:7]           
    lattice = np.float64([a_vec, b_vec, c_vec])
    atom_type, atom_position = get_atoms(outfile)
    atom_number = cell_utils.convert_atomtype(atom_type)
    cell = (lattice, atom_position, atom_number)
    
    out = {}
    out['atom'] = atom_type
    out['cell'] = cell
    out['efermi'] = get_value(outfile, keyword="Fermi energy")
    out['kpts'] = 0 # kmesh info: TODO: will test the Gamma point case first

    return out     
        
def read_bs(filename):
    '''Read the cp2k *.bs file'''
    assert check_exist(filename), "Cannot find " + filename

    with open(filename, "r") as bs_file:
        symm_k_coor = []
        kpath_frac = []
        bands = []
        for kpoint_set in SET_MATCH.finditer(bs_file.read()):
            nkpts = np.int64(kpoint_set.groupdict()['totalpoints'])
            set_content = kpoint_set.group('content')
            sym_kpoint_coor = []
            for point in SPOINTS_MATCH.finditer(set_content):
                point_dict = point.groupdict()
                sym_kpoint_coor.append([point_dict['a'], point_dict['b'], point_dict['c']])
                
            symm_k_coor.append(np.float64(sym_kpoint_coor))
            kpts = []
            band_alpha = []
            band_beta = []
            for point in POINTS_MATCH.finditer(set_content):
                point_dict = point.groupdict()
                spin = point_dict['spin']
                if spin == '1':
                    spin_polarized = False
                    energies = np.float64(point_dict['values'].split())
                    band_alpha.append(energies)
                    kpt = [point_dict['a'], point_dict['b'], point_dict['c']]
                    kpts.append(kpt)
                else:
                    spin_polarized = True
                    energies = np.float64(point_dict['values'].split())
                    band_beta.append(energies)  
            band_alpha
            kpath_frac.append(np.float64(kpts))
            band = [band_alpha]
            if spin_polarized is True: 
                band.append(band_beta)   
            bands.append(np.asarray(band))
    
    out = {}
    out['band'] = bands
    out['symm_k_coor'] = symm_k_coor
    out['kpath_frac'] = kpath_frac 

    return out     
        

        