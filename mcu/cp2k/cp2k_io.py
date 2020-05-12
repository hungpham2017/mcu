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


# Get no. of electron
def read_cell(data, key):
    '''CP2K has different CELL block with similar template
       This function get cell starting with the provided key, such as: CELL, CELL_TOP, CELL_REF, CELL_UC'''
    CELL_MATCH = re.compile(r'''
    [\w\W]* 
    ''' + key + '''\| [ ]* Vector [ ]* a [ ]* [angstrom [\]:]* (?P<a1>\S+) [ ]* (?P<a2>\S+) [ ]* (?P<a3>\S+) [ ]* \|a\| [ ]* = [ ]* (?P<a>\S+) [ \n]* 
    ''' + key + '''\| [ ]* Vector [ ]* b [ ]* [angstrom [\]:]* (?P<b1>\S+) [ ]* (?P<b2>\S+) [ ]* (?P<b3>\S+) [ ]* \|b\| [ ]* = [ ]* (?P<b>\S+) [ \n]*
    ''' + key + '''\| [ ]* Vector [ ]* c [ ]* [angstrom [\]:]* (?P<c1>\S+) [ ]* (?P<c2>\S+) [ ]* (?P<c3>\S+) [ ]* \|c\| [ ]* = [ ]* (?P<c>\S+)   
    ''', re.VERBOSE)
    
    cell_data = CELL_MATCH.match(data)
    if cell_data is not None:
        a = [cell_data['a1'], cell_data['a2'], cell_data['a3']]
        b = [cell_data['b1'], cell_data['b2'], cell_data['b3']]
        c = [cell_data['c1'], cell_data['c2'], cell_data['c3']]
        return np.float64([a, b, c])
    else:
        return None
        
atoms_MATCH = re.compile(r'''
[\w\W]*
  Atom [ ]* Kind [ ]* Element [ ]* X [ ]* Y [ ]* Z [ ]* Z\(eff\) [ ]* Mass
  \n\n
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)
         
atom_MATCH = re.compile(r'''
[ ]*
 (?P<order>\d+) [ ]* (?P<kind>\d+) [ ]* (?P<element>\S+) [ ]*  (?P<Z>\d+) [ ]* (?P<x>\S+) [ ]* (?P<y>\S+) [ ]* (?P<z>\S+) [ ]* (?P<Zeff>\S+) [ ]* (?P<mass>\S+)
''', re.VERBOSE)

def read_out(filename):  
    '''Read the cp2k prefix.out file'''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, "r") as data_file:
        data = data_file.read()

        # Cell information        
        lattice = read_cell(data, 'CELL')
        
        # Get atom info
        atoms_blocks = []
        for atoms in atoms_MATCH.finditer(data):
            atoms_data = atoms.groupdict()['content']
            atom = []
            abs_coords = []
            Z = []
            Zeff = []
            mass = []
            for atm in atom_MATCH.finditer(atoms_data):
                atom_dict = atm.groupdict()
                Z.append(atom_dict['Z'])
                Zeff.append(float(atom_dict['Zeff']))
                mass.append(atom_dict['mass'])
                atom.append(atom_dict['element'])
                abs_coords.append([atom_dict['x'], atom_dict['y'], atom_dict['z']])
            
            abs_coords = np.float64(abs_coords)
            frac_coords = abs_coords @ np.linalg.inv(lattice)
            
            atoms = {}
            atoms['Z'] = Z
            atoms['Zeff'] = Zeff
            atoms['mass'] = mass
            atoms['atom'] = atom
            atoms['abs_coords'] = abs_coords
            atoms['frac_coords'] = frac_coords
            atoms_blocks.append(atoms)
         
        # Taking info from the last geometry block 
        atoms = atoms_blocks[-1]
        atom_number = cell_utils.convert_atomtype(atoms['atom'])
        cell = (lattice, atoms['frac_coords'], atom_number)
        
        out = {}
        out['atoms blocks'] = atoms_blocks
        out['atom'] = atoms['atom']
        out['Zeff'] = atoms['Zeff']
        out['cell'] = cell
        out['kpts'] = 0 # kmesh info: TODO: will test the Gamma point case first

        return out             
    
#--------------------------------------------------------------------------------------------------------- 
# the read_bs function was modified from cp2k_bs2csv.py
# ref: https://www.cp2k.org/exercises:2016_uzh_cmest:band_structure_calculation
#---------------------------------------------------------------------------------------------------------
         
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
  POINT [ ]+ (?P<pointnr>\d+) [ ]+ (?P<klabel>\S*) [ ]+ (?P<a>\S+) [ ]+ (?P<b>\S+) [ ]+ (?P<c>\S+)
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
        
def read_bs(filename):
    '''Read the cp2k *.bs file'''
    assert check_exist(filename), "Cannot find " + filename

    with open(filename, "r") as data_file:
        data = data_file.read()
        symm_k_coor = []
        symm_k_label = []
        kpath_frac = []
        bands = []
        for kpoint_set in SET_MATCH.finditer(data):
            nkpts = np.int64(kpoint_set.groupdict()['totalpoints'])
            set_content = kpoint_set.group('content')
            sym_kpoint_coor = []
            sym_kpoint_label = []
            for point in SPOINTS_MATCH.finditer(set_content):
                point_dict = point.groupdict()
                sym_kpoint_label.append(point_dict['klabel'])
                sym_kpoint_coor.append([point_dict['a'], point_dict['b'], point_dict['c']])
                
            symm_k_coor.append(np.float64(sym_kpoint_coor))
            symm_k_label.append(sym_kpoint_label)
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
    out['symm_k_label'] = symm_k_label
    out['kpath_frac'] = kpath_frac 

    return out     
        

        