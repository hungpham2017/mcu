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
from mcu.utils.misc import check_exist
from mcu.cell import utils  

#--------------------------------------------------------------------------------------------------------- 
# USEFUL FUNCTIONS TO MANIPULATE wannier90 files
#---------------------------------------------------------------------------------------------------------
def copy_block(blocks, keyword, convert_to_float=False):
    ''' copy info in a block'''
    start_key = ('begin ' + keyword).lower()
    end_key = ('end '+ keyword).lower()  
    block = []
    copy = False
    for line in blocks:
        if end_key == line[:len(end_key)].lower():            
            copy = False   
            
        if copy == True: 
            line = line.replace('\n','')
            if convert_to_float == True:
                line = np.float64(line.split())
            block.append(line)
        else:
            if start_key == line[:len(start_key)].lower(): copy = True
            
    return block 
    
def extract_parameter(blocks, keyword):
    '''Get the parameter after a keyword'''
    value = None
    for line in blocks:
        if keyword.lower() in line.lower():  
            return line.split()[2:]
            

    
#--------------------------------------------------------------------------------------------------------- 
# MAIN CLASS
#---------------------------------------------------------------------------------------------------------
class io:
    def __init__(self, seedname="wannier90"):
        '''
            seedname can contain the path as well if the wannier90 files are not in the working directory
        '''
        self.seedname = seedname

    def read_win(self, seedname=None):  
        '''Read the seedname.win file'''
        
        if seedname is None: seedname = self.seedname
        if check_exist(seedname + '.win'):
            win = open(seedname + '.win', "r").readlines() 

            # Cell information
            lattice = copy_block(win, 'Unit_Cell_Cart', convert_to_float=True)
            atom_block = copy_block(win, 'atoms_cart')
            atom_sym = []
            atom_position = []
            for atm in atom_block:
                temp = atm.split()
                atom_sym.append(temp[0])
                atom_position.append(np.float64(temp[1:]))
            self.atom = atom_sym
            atom_number = utils.convert_atomtype(self.atom )
            self.cell = (lattice, atom_position, atom_number)
            
            # kmesh info
            self.kpts = np.asarray(copy_block(win, 'kpoints', convert_to_float=True))
            self.klabel = None
            kpoint_path = copy_block(win, 'kpoint_path')
            if kpoint_path is not []:
                num_kpoint = len(kpoint_path)
                high_sym_kpoints = []
                for i, kpoint in enumerate(kpoint_path):
                    if i == (num_kpoint - 1):
                        temp = kpoint.split()
                        high_sym_kpoints.append([temp[0], np.float64(temp[1:4])])
                        high_sym_kpoints.append([temp[4], np.float64(temp[5:])])
                    else:
                        temp = kpoint.split()
                        high_sym_kpoints.append([temp[0], np.float64(temp[1:4])])
                self.klabel = high_sym_kpoints
            
        else:
            print('Cannot find the *.win file. Check the path:', seedname + '.win')
        
    def read_wout(self, seedname=None):  
        '''Read the seedname.win file'''
        if seedname is None: seedname = self.seedname
        
    def read_band(self, seedname=None):  
        '''Read the seedname_band.dat file'''
        if seedname is None: seedname = self.seedname
        if check_exist(seedname + '_band.dat'):
            with open(seedname + '_band.dat', "r") as file:
                lines = file.read().splitlines()
            
            band = []
            for i, line in enumerate(lines):
                if line == '  ':
                    nkpts = i + 1
                    break
                else:
                    band.append(np.float64(line.split()))
                    
            nbands = len(lines) // nkpts
            bands = [band]
            for line in range(1,nbands):
                band = []
                for point in range(nkpts):
                    band.append(lines[line*nkpts + point].split())
                bands.append(band[:-1])
            bands = np.asarray(bands, dtype=np.float64)
            self.kpath = bands[0][:,0]
            self.band = bands[:,:,1].T
                
        else:
            print('Cannot find the *_band.dat file. Check the path:', seedname + '_band.dat')
        



        