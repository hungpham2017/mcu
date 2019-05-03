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
import subprocess
from mcu.utils import misc
from mcu.cell import utils
        
class cif:
    def __init__(self, file=None):
    
        if file == None:        # The 1st 
            proc = subprocess.Popen('/bin/ls *.cif', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = proc.communicate()
            file = str(out).split("'")[1].split("\\")[0]
        if not misc.check_exist(file):
            print('Cannot find any cif file in the current directory', file)
        else:
            self.cif = open(file, "r").readlines()
            
    def read_lattice(self):
        '''Read lattice constant'''
        
        a = self.extract_line('_cell_length_a', data_type='float')
        b = self.extract_line('_cell_length_b', data_type='float')
        c = self.extract_line('_cell_length_c', data_type='float')
        alpha = self.extract_line('_cell_angle_alpha', data_type='float')
        beta = self.extract_line('_cell_angle_beta', data_type='float')
        gamma = self.extract_line('_cell_angle_gamma', data_type='float')
        sg_sym = self.extract_line('_symmetry_space_group_name_H-M', data_type='str')
        sg_num = self.extract_line('_symmetry_Int_Tables_number', data_type='int') 
        self.symbol, self.frac = self.extract_coordinate()
        
        
    def extract_coordinate(self):
        '''Extract the (irreducible) coordinates block'''
        out = None
        block_keys = []
        start = 0
        stop = False
        for i in range(len(self.cif)):
            line = self.cif[i]
            if line.startswith('_atom_site'):
                block_keys.append(line.split('_atom_site_')[1].split('\n')[0])
                start = i + 1
                stop = True
            elif stop == True:
                break
                
        num_keys = len(block_keys)

        where_symbol = block_keys.index('type_symbol')
        where_x = block_keys.index('fract_x')
        where_y = block_keys.index('fract_y')
        where_z = block_keys.index('fract_z')
        if 'label' in block_keys: 
            where_label = block_keys.index('label')
            label = []
        if 'occupancy' in block_keys:
            where_occ = block_keys.index('occupancy')
            occupancy = []
            
        symbol = []
        frac = []
        for i in range(start, len(self.cif)):
            line = self.cif[i].split()
            if len(line) == num_keys:
                symbol.append(line[where_symbol])
                frac.append([utils.rm_paren(line[where_x]),
                             utils.rm_paren(line[where_y]),
                             utils.rm_paren(line[where_z]]))
                if 'label' in block_keys: label.append(line[where_label])
                if 'occupancy' in block_keys: occupancy.append(line[where_occ])                
            else:
                break

        frac = np.float64(frac)
        
        return symbol, frac
        
    def extract_sym_operator(self):
        '''Extract the (irreducible) coordinates block'''
        
        start = 0
        for i in range(len(self.cif)):
            line = self.cif[i]
            if line.startswith('_symmetry_equiv_pos_as_xyz'):
                start = i + 1
                break
                
        sym_operators = []
        for i in range(start, len(self.cif)):
            line = self.cif[i]
            if line.strip() == '' or line.startswith("_"): 
                break
                    
            if "'" in line:
                line = line.split("'")[1]
            else:
                line = line.split()[0]
               
            sym_operators.append(line)
            
        return sym_operators
                
                
    def extract_line(self, key, data_type='str'):
        '''Extract the value after a keyword'''
        out = None
        for line in self.cif:
            if line.startswith(key):
                temp = line.split()
                if len(temp) == 2:
                    out = utils.rm_paren(temp[1].replace("'",""))
                else:
                    new_temp = ''
                    for i in range(len(temp)-1):
                        new_temp += temp[i+1] + ' '
                    out = new_temp.replace("'","").strip()
                break
                
        if out != None:
            if data_type == 'float':
                out = float(out)
            elif data_type == 'int':
                out = int(out)
            
        return out
            
            
                

               
        
        
                
                
                
                
            

    
        
        

        

        