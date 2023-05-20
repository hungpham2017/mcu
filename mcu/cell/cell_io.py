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
from ..utils import misc
from . import utils as cell_utils
from . import spg_wrapper


'''=======================================       
 I/O for CIF
=======================================''' 

def get_value(data, key):
    '''Giving a keyword, return the value next to it'''
    key = key.strip()
    keyword_MATCH = re.compile(r'''[\w\W]* ''' + key + ''' [ ]* (?P<value>\S+)''', re.VERBOSE)
    match = keyword_MATCH.match(data)
    if match is None:
        return match
    else:
        return match['value']


symmetry1_loop_MATCH = re.compile(r'''
[\w\W]* loop_  [ ]* [\n]* [ ]* _symmetry_equiv_pos_as_xyz 
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|_|loop_)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

symmetry2_loop_MATCH = re.compile(r'''
[\w\W]* loop_  [ ]* [\n]* [ ]* _symmetry_equiv_pos_site_id [\n]* _symmetry_equiv_pos_as_xyz 
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|_|loop_)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

def find_atom_site_string(data):
    temp = data.split()
    nitems = len(temp)
    if '_atom_site_fract_x' in temp:
        fract_x_idx = temp.index('_atom_site_fract_x')
        left_list = []
        for i in range(fract_x_idx, 0, -1):
            if '_atom_site_' in temp[i]:
                left_list.append(temp[i])
            else:
                break
        right_list = []
        for i in range(fract_x_idx, nitems, 1):
            if '_atom_site_' in temp[i]:
                right_list.append(temp[i])
            else:
                break
        left_list.reverse()
        return  left_list[:-1] + right_list
    else:
        None
           
def read_cif(filename):  
    '''Read cif file'''

    assert misc.check_exist(filename), 'Cannot find : ' + filename
    with open(filename, 'r') as data_file:
        data = data_file.read()
        
        # get lattice constans
        a = get_value(data, '_cell_length_a').replace("(","").replace(")","")
        b = get_value(data, '_cell_length_b').replace("(","").replace(")","")
        c = get_value(data, '_cell_length_c').replace("(","").replace(")","")     
        alpha = get_value(data, '_cell_angle_alpha').replace("(","").replace(")","")
        beta = get_value(data, '_cell_angle_beta').replace("(","").replace(")","")        
        gamma = get_value(data, '_cell_angle_gamma').replace("(","").replace(")","")        
        lattice = np.float64([a, b, c, alpha, beta, gamma])
        
        
        # Get symmetry
        space_group_name = get_value(data, '_symmetry_space_group_name_H-M')
        if space_group_name is not None: 
            space_group_name = space_group_name.replace("'","")
        if space_group_name is not None: 
            space_group_name = space_group_name.strip().replace("'","")  
            if space_group_name == "?": space_group_name = None
            
        space_group_number = get_value(data, '_symmetry_Int_Tables_number')
        cell_setting = get_value(data, '_symmetry_cell_setting')
        
        symmetry_operatator = None
        if symmetry1_loop_MATCH.match(data) is not None:
            symmetry_loop = symmetry1_loop_MATCH.match(data)
            if "'" in symmetry_loop['content']:
                temp = re.findall(r"""'(.*?)'""", symmetry_loop['content'])
                symmetry_operatator = [item.replace(" ","") for item in temp]
            else:
                symmetry_operatator = symmetry_loop['content'].split()
            
        elif symmetry2_loop_MATCH.match(data) is not None:
            symmetry_loop = symmetry2_loop_MATCH.match(data)
            if "'" in symmetry_loop['content']:
                temp = re.findall(r"""'(.*?)'""", symmetry_loop['content'])
                symmetry_operatator = [item.replace(" ","") for item in temp]
            else:
                temp = symmetry_loop['content'].split()
                symmetry_operatator = [temp[2*item + 1] for item in range(len(temp)//2)]
         
        # Get atom positions and labels
        
        found_atoms = False
        atom_site_keys = find_atom_site_string(data)

        atom_loop = None
        if atom_site_keys is not None:
            atom_site_string = "[ ]* [\n]* [ ]*".join(atom_site_keys)
            atom_loop_MATCH = re.compile(r'''
            [\w\W]* loop_ [ ]* [\n]* [ ]* ''' + atom_site_string +''' 
            (?P<content>
              [\s\S]*?(?=\n.*?[ ] $|$|loop_) 
            )
            ''', re.VERBOSE)
            atom_loop = atom_loop_MATCH.match(data)

        fract_xyz = None
        if atom_loop is not None:
            found_atoms = True
            x_idx = atom_site_keys.index('_atom_site_fract_x')  
            y_idx = atom_site_keys.index('_atom_site_fract_y')
            z_idx = atom_site_keys.index('_atom_site_fract_z')
            
            # Check all the atom_site keyword
            if '_atom_site_label' in atom_site_keys: 
                label_idx = atom_site_keys.index('_atom_site_label')
            else:
                label_idx = None
                
            if '_atom_site_type_symbol' in atom_site_keys: 
                symbol_idx = atom_site_keys.index('_atom_site_type_symbol') 
            else:
                symbol_idx = None
                
            if '_atom_site_occupancy' in atom_site_keys: 
                occupancy_idx = atom_site_keys.index('_atom_site_occupancy') 
            else:
                occupancy_idx = None

            symbols = []
            labels = []
            fract_xyz = []
            occupancy = []
            atom_content = atom_loop['content']
            atom_string = "".join(['[\n]* [ ]* (?P<value%d>\S+) [ ]* '% (i) for i in np.arange(len(atom_site_keys))]) 
            atom_MATCH = re.compile(atom_string, re.VERBOSE)
            for i, atom in enumerate(atom_MATCH.finditer(atom_content)):
                atom_dict = atom.groupdict()

                x_string = atom_dict['value' + str(x_idx)].replace("(","").replace(")","")
                y_string = atom_dict['value' + str(y_idx)].replace("(","").replace(")","")
                z_string = atom_dict['value' + str(z_idx)].replace("(","").replace(")","")
                
                def is_float(string):
                    return string.replace(".","").replace("-","").isdigit()
                    
                if is_float(x_string) and is_float(y_string) and is_float(z_string):
                    x = float(x_string)
                    y = float(y_string)
                    z = float(z_string)
                else:    
                    break
                
                fract_xyz.append([x, y, z])
                
                if symbol_idx is None:
                    symbols.append(None)
                else:    
                    symbols.append(atom_dict['value' + str(symbol_idx)].capitalize())
                
                if label_idx is None:
                    labels.append(None)
                else:    
                    labels.append(atom_dict['value' + str(label_idx)])
                    
                if occupancy_idx is None:
                    occupancy.append(1.0)
                else:    
                    occ = atom_dict['value' + str(occupancy_idx)].replace("(","").replace(")","")
                    occupancy.append(float(occ))
                    
            if occupancy_idx is not None:
                for i, occ in enumerate(occupancy):
                    if abs(occ - 1.0) > 1.e-4:
                        print("WARNING! {0:6s} {1:4s} has a fractional occupancy of {2:2.4f}".format(str(labels[i]), str(symbols[i]), occ))
                    
            fract_xyz = np.float64(fract_xyz)
        
        out = {}
        out['lattice parameters'] = lattice
        out['spg number'] = space_group_number
        out['spg name'] = space_group_name
        out['symmetry operators'] = symmetry_operatator
        out['cell_setting'] = cell_setting
        if found_atoms:
            out['fract_xyz'] = fract_xyz
            out['symbols'] = symbols    # Must be an elemenent name
            out['labels'] = labels      # Can be any string
            out['occupancy'] = occupancy 
        else:
            out['fract_xyz'] = None
            out['symbols'] = None       # Must be an elemenent name
            out['labels'] = None        # Can be any string
            out['occupancy'] = None
        return out


def cif2cell(filename):
    '''Read a cif file and return a cell spglib object'''
    assert misc.check_exist(filename), "Cannot find : " + filename
    
    out = read_cif(filename)
    lattice = cell_utils.convert_lattice(out['lattice parameters'])
    assert out['fract_xyz'] is not None, "The provided CIF does not have atom information. Check your CIF: " + filename
    irred_symbol = out['symbols']
    irred_frac = out['fract_xyz']
    irred_occupancy = out['occupancy'] 
    if int(sum(irred_occupancy)) != len(irred_occupancy):
        print("WARNING! Cannot simply convert a disordered CIF to a cell object. Refine your CIF")
        return None
    else:
        symopt = out['symmetry operators']
        rotations, translations = cell_utils.symop_xyz2mat(symopt)
        full_symbol, full_frac = cell_utils.atoms_irred2full(lattice, irred_symbol, irred_frac, rotations, translations)
        if full_symbol is None:
            print("WARNING! Cannot simply convert a disordered CIF to a cell object. Refine your CIF")
            return None
        else:
            atom_type = cell_utils.convert_atomtype(full_symbol)
            return (lattice, full_frac, atom_type)

 
'''=======================================       
 EXPORT: CIF, XSF, POSCAR 
======================================='''         
def write_poscar(cell, filename=None):
    comment = misc.date()
    lattice = np.asarray(cell[0])
    positions = np.asarray(cell[1])
    atoms = np.asarray(cell[2])
    idx = np.argsort(atoms)
    atoms = atoms[idx]
    symbol = cell_utils.convert_atomtype(atoms)
    irred_symbol, count = misc.unique(symbol) 
    positions = positions[idx]
    
    with open(filename, 'w') as f:
        f.write('Generated by mcu on ' + comment + '\n')
        f.write('1.0\n') 
        for i in range(3):
            f.write('   %15.10f %15.10f %15.10f\n' % (lattice[i,0],lattice[i,1],lattice[i,2]))         
        for symb in irred_symbol:
            f.write(' %s ' % (symb)) 
        f.write('\n')
        for num_atom in count:
            f.write(' %d ' % (num_atom)) 
        f.write('\n')
        f.write('Direct\n')
        for atom in positions:
            f.write('   %15.10f %15.10f %15.10f\n' % (atom[0],atom[1],atom[2]))   

def write_xsf(cell, filename=None):
    comment = misc.date()
    lattice = np.asarray(cell[0])
    positions = np.asarray(cell[1])
    abs_positions = positions.dot(lattice)
    atoms = np.asarray(cell[2])
    natom = len(atoms)
    symbol = cell_utils.convert_atomtype(atoms)

    with open(filename + '.xsf', 'w') as f:
        f.write('Generated by mcu on ' + comment + '\n')   
        f.write('CRYSTAL\n')
        f.write('PRIMVEC\n')    
        for i in range(3):
            f.write('   %15.10f %15.10f %15.10f\n' % (lattice[i,0],lattice[i,1],lattice[i,2]))   
        f.write('CONVVEC\n')
        for i in range(3):
            f.write('   %15.10f %15.10f %15.10f\n' % (lattice[i,0],lattice[i,1],lattice[i,2]))     
        f.write('PRIMCOORD\n')
        f.write('%3d %3d\n' % (natom, 1))
        for atom in range(natom):
            f.write(' %s  %15.10f  %15.10f  %15.10f\n' % (symbol[atom], abs_positions[atom][0], abs_positions[atom][1], abs_positions[atom][2]))          
      
def write_cif(cell, space_group, symopt, filename=None):
    if filename == None: filename = 'structure_mcu'
    comment = misc.date()
    lattice_mat = np.asarray(cell[0])
    lattice = cell_utils.convert_lattice(lattice_mat)   
    positions = np.asarray(cell[1])
    atoms = np.asarray(cell[2])
    natom = len(atoms)
    symbol = cell_utils.convert_atomtype(atoms)
    nsymopt = len(symopt)
    
    with open(filename + '.cif', 'w') as f:
        f.write('data_New_Crystal\n') 
        f.write("_audit_creation_method         '%s'\n" % ('Generated by mcu on ' + comment))
        f.write('_cell_length_a     %15.10f\n' % (lattice[0]))
        f.write('_cell_length_b     %15.10f\n' % (lattice[1])) 
        f.write('_cell_length_c     %15.10f\n' % (lattice[2]))   
        f.write('_cell_angle_alpha     %15.10f\n' % (lattice[3]))
        f.write('_cell_angle_beta      %15.10f\n' % (lattice[4])) 
        f.write('_cell_angle_gamma     %15.10f\n' % (lattice[5])) 
        f.write('\n')  
        f.write("_symmetry_space_group_name_H-M     '%s'\n" % (space_group[1]))
        f.write('_symmetry_Int_Tables_number        %d\n' % (space_group[0]))
        f.write('loop_\n')        
        f.write('_symmetry_equiv_pos_as_xyz\n') 
        for i in range(nsymopt):
            f.write('%s\n' % (symopt[i]))
            
        f.write('\n')  
        f.write('loop_\n')
        f.write('_atom_site_label\n')
        f.write('_atom_site_type_symbol\n')
        f.write('_atom_site_fract_x\n')
        f.write('_atom_site_fract_y\n')
        f.write('_atom_site_fract_z\n')        
        for atom in range(natom):
            f.write('   %s   %s   %15.10f   %15.10f   %15.10f\n' % (symbol[atom], symbol[atom], positions[atom][0], positions[atom][1], positions[atom][2]))     

            

