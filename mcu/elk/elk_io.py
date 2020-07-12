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
from ..vasp import const
from ..cell import utils as cell_utils  


lattice_MATCH = re.compile(r'''
[\w\W]*
 Lattice [ ]* vectors [ ]* \:
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n[ ]\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)
         
species_MATCH = re.compile(r'''
[ ]*
 Species [ ]* \:
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n[ ]\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

kmesh_MATCH = re.compile(r'''
[\w\W]* k-point [ ]* grid [ ]* \: [ ]* (?P<kx>\S+) [ ]* (?P<ky>\S+) [ ]* (?P<kz>\S+)
''', re.VERBOSE)

koffset_MATCH = re.compile(r'''
[\w\W]* k-point [ ]* offset [ ]* \: [ ]* (?P<kx>\S+) [ ]* (?P<ky>\S+) [ ]* (?P<kz>\S+)
''', re.VERBOSE)

total_charge_MATCH = re.compile(r'''
[\w\W]*
Total [ ]* nuclear [ ]* charge    [ ]* \: [ ]* (?P<nuclear>\S+) [ ]* \n  
Total [ ]* core [ ]* charge       [ ]* \: [ ]* (?P<core>\S+) [ ]* \n  
Total [ ]* valence [ ]* charge    [ ]* \: [ ]* (?P<valence>\S+) [ ]* \n     
Total [ ]* excess [ ]* charge     [ ]* \: [ ]* (?P<excess>\S+) [ ]* \n     
Total [ ]* electronic [ ]* charge [ ]* \: [ ]* (?P<electronic>\S+)   
''', re.VERBOSE)

states_charge_MATCH = re.compile(r'''
[\w\W]*
Number [ ]* of [ ]* empty [ ]* states  [ ]* \: [ ]* (?P<empty>\S+) [ ]* \n
Total [ ]* number [ ]* of [ ]* valence [ ]* states [ ]* \: [ ]* (?P<valence>\S+) [ ]* \n
Total [ ]* number [ ]* of [ ]* core [ ]* states [ ]* \: [ ]* (?P<core>\S+)
''', re.VERBOSE)

############ For SCF loops ################## 
energies_MATCH = re.compile(r'''
[ ]*
 Energies [ ]* \:
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n[ ]\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

charges_MATCH = re.compile(r'''
[ ]*
 Charges [ ]* \:
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n[ ]\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

def read_info(filename):  
    '''Read the INFO.OUT file'''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, "r") as data_file:
        data = data_file.read()

        # spin
        spin = int('spin-polarised' in data)
        soc = 'spin-orbit coupling' in data
        if soc: spin = False

        # Cell information
        # row-vectors instead of column-vectors as in elk output
        tmp = lattice_MATCH.match(data)['content']
        lattice = np.float64(tmp.split()).reshape(3,3, order='F')      
        
        # Get atom info
        species_blocks = []
        atom = []
        magnetic = []
        frac_coords = []
        for atoms in species_MATCH.finditer(data):
            atoms_data = atoms.groupdict()['content']
            tmp = atoms_data.split('\n')
            species = {}
            species['atom type'] = re.findall(r"\((?P<val>\S+)\)", tmp[0])[0]
            species['parameters'] = tmp[1].split(':')[1].strip()
            species['name'] = tmp[2].split(':')[1].strip()
            species['nuclear charge'] = float(tmp[3].split(':')[1])
            species['electronic charge'] = float(tmp[4].split(':')[1])
            species['atomic mass'] = float(tmp[5].split(':')[1])
            species['muffin-tin radius'] = float(tmp[6].split(':')[1])
            species['muffin-tin radial'] = int(tmp[7].split(':')[1])
            species['muffin-tin inner'] = int(tmp[8].split(':')[1])
            for atm in tmp[10:]:
                atom.append(species['atom type'])
                atm_split = atm.split()
                frac_coords.append(np.float64(atm_split[2:5]))
                magnetic.append(np.float64(atm_split[5:]))

            species_blocks.append(species)
         
        # Make the cell tuple
        atom_number = cell_utils.convert_atomtype(atom)
        cell = (lattice * const.AUTOA, frac_coords, atom_number)
        
        # kmesh
        tmp = kmesh_MATCH.match(data)
        kmesh = np.int64([tmp['kx'], tmp['ky'], tmp['kz']])
        tmp = koffset_MATCH.match(data)
        kmesh_shift = np.float64([tmp['kx'], tmp['ky'], tmp['kz']])
        
        # total charge
        charge_data = total_charge_MATCH.match(data)
        charge = {}
        charge['nuclear'] = float(charge_data['nuclear'])
        charge['core'] = float(charge_data['core'])
        charge['valence'] = float(charge_data['valence'])
        charge['excess'] = float(charge_data['excess'])
        charge['electronic'] = float(charge_data['electronic'])
        
        # Total states:
        states_data = states_charge_MATCH.match(data)
        states = {}
        states['empty'] = states_data['empty']
        states['core'] = states_data['core']
        states['valence'] = states_data['valence']
        
        
        # Collect the SCF cycles
        energies_blocks = []
        for cycle in energies_MATCH.finditer(data):
            energies_data = cycle.groupdict()['content']
            tmp = energies_data.split('\n')
            energies = {}
            energies['Fermi'] = float(tmp[1].split(':')[1])
            energies['sum of eigenvalues'] = float(tmp[2].split(':')[1])
            energies['electron kinetic'] = float(tmp[3].split(':')[1])
            energies['core electron kinetic'] = float(tmp[4].split(':')[1])
            energies['Coulomb'] = float(tmp[5].split(':')[1])
            energies['Coulomb potential'] = float(tmp[6].split(':')[1])
            energies['nuclear-nuclear'] = float(tmp[7].split(':')[1])
            energies['electron-nuclear'] = float(tmp[8].split(':')[1])
            energies['Hartree'] = float(tmp[9].split(':')[1])
            energies['Madelung'] = float(tmp[10].split(':')[1])
            energies['xc potential'] = float(tmp[11].split(':')[1])
            energies['exchange'] = float(tmp[12].split(':')[1])
            energies['correlation'] = float(tmp[13].split(':')[1])
            energies['electron entropic'] = float(tmp[14].split(':')[1])
            energies['total energy'] = float(tmp[15].split(':')[1])
            energies_blocks.append(energies)        

        out = {}
        out['Species'] = species_blocks
        out['Energies'] = energies_blocks
        out['charge'] = charge
        out['states'] = states
        out['atom'] = atom
        out['cell'] = cell
        out['soc'] = soc
        out['spin'] = spin
        
        return out             

def read_kpoints(filename):  
    '''Read the KPOINTS.OUT file'''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, "r") as data_file:
        data = data_file.read()
        tmp = np.float64(" ".join(data.split('\n')[1:-1]).split()).reshape(-1, 6)
        return tmp[:,1:4], tmp[:,4]           # kpts and kpts weight
        
eigval_MATCH = re.compile(r'''
[ ]*
 \([ ]* state, [ ]* eigenvalue [ ]* and [ ]* occupancy [ ]* below\)
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n[ ]\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

def read_eigval(filename):  
    '''Read the KPOINTS.OUT file
       Return:
        eigvals[spin, kpt, band, eigenvalues]
        occ[spin, kpt, band, occupancy]
    '''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, "r") as data_file:
        data = data_file.read()

        eigvals = []
        occ = []
        for eigval in eigval_MATCH.finditer(data):
            eigval_data = eigval['content']
            tmp = np.float64(eigval_data.split()).reshape(-1, 3)
            eigvals.append(tmp[:,1])
            occ.append(tmp[:,2])
            
        eigvals = np.asarray(eigvals)
        nkpts = eigvals.shape[0]
        nband = eigvals.shape[1]
        occ = np.asarray(occ)
        
        # Check if this is a spin-polarized calculation 
        if (occ[0] <= 1.0).any() and abs((-occ[0]).argsort() - np.arange(occ[0].shape[0])).sum() > 1e-3: 
            assert np.mod(nband,2) == 0, "This is not a spin-polarised calculation"
            nband = nband//2
        
        eigvals = eigvals.reshape(nkpts, nband, -1).transpose(2, 0, 1) 
        occ = occ.reshape(nkpts, nband, -1).transpose(2, 0, 1) 
        return eigvals, occ 
        
############ Read BAND***.OUT - begin ##################        
def read_bandlines(filename):  
    '''Read the BANDLINES.OUT file
    '''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, "r") as data_file:
        data = np.float64(data_file.read().split()).reshape(-1,2)
        nkpoint = data.shape[0] // 2
        return data[2*np.arange(nkpoint), 0] / const.AUTOA

def read_band(filename, spin_polarized=False):
    '''Read BAND.OUT or BAND_Sss_Aaaaa.OUT

       Depends on the tasks, BAND.OUT or BAND_Sss_Aaaaa.OUT contains different projection
       For tasks 20:    kpt  |  e (a.u.)
       For tasks 21:    kpt  |  e (a.u.)  |  Total  |   s   |   p   |   d   |   f   |
       For tasks 22:    kpt  |  e (a.u.)  |   s     |  p(-1,0,+1)  |  d(-2,-1,0,+1,+2)  |  f(-3,-2,-1,0,+1,+2,+3)
       For tasks 23:    kpt  |  e (a.u.)  |   up    |  down |
       
       Return:
                bands[spin, kpt, band-th, projected contribution] 
    '''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, 'r') as data_file:
        data = data_file.read()
        tmp = data.split('\n \n')[:-1]
        
        # Figure out how many kpt and how many column
        nband = len(tmp)
        tmp2 = tmp[0].split('\n')
        nkpts = len(tmp2)
        ncol = len(tmp2[0].split())
        bands = []
        for i, band in enumerate(tmp):
            formatted_band = np.float64(band.split()).reshape(-1, ncol)    
            if i == 0: proj_kpath = formatted_band[:,0] / const.AUTOA           # in Angstrom-1
            bands.append(formatted_band[:,1:]) 

        if spin_polarized: 
            assert np.mod(nband,2) == 0, "This is not a spin-polarised calculation"
            nband = nband//2
            
        bands = np.asarray(bands).reshape(-1, nband, nkpts, ncol - 1).transpose(0, 2, 1, 3) * const.AUTOEV   
        return proj_kpath, bands 

def read_plot1d(filename):
    '''Read the plot1d of elk.in'''
    
    with open(filename, 'r') as data_file:
        data = data_file.read()
        assert 'plot1d' in data, "plot1d block cannot be found in the elk.in" 
        tmp_with_comments = data.split('\n')
        tmp = [item for item in tmp_with_comments if not '!' in item]   # remove comments
        for i, item in enumerate(tmp):
            if 'plot1d' in item:
                where_plot1d = i
                break
     
        nkpoint, npts = np.int64(tmp[where_plot1d + 1].split()[:2])
        sym_kvecs = np.float64([tmp[i].split()[:3] for i in range(where_plot1d + 2, where_plot1d + nkpoint + 2)])
        
        return sym_kvecs, npts
        
############ Read DOS - begin ##################  
def read_dos(filename, spin_polarized=False):
    '''Read TDOS.OUT, IDOS.OUT, PDOS_Sss_Aaaaa.OUT

       For PDOS_Sss_Aaaaa.OUT:       
       e (a.u.)  |   s     |  p(-1,0,+1)  |  d(-2,-1,0,+1,+2)  |  f(-3,-2,-1,0,+1,+2,+3)
       Return:
                dos[spin, e, dos] 
    '''
    assert check_exist(filename), 'Cannot find : ' + filename
    with open(filename, 'r') as data_file:
        data = data_file.read()
        tmp = data.split('\n \n')[:-1]
        
        # Figure out how many kpt and how many column
        ndos = len(tmp)
        nepsilon = len(tmp[0].split('\n'))
        dos = []
        for i, block in enumerate(tmp):
            formatted_dos = np.float64(block.split()).reshape(nepsilon, 2)    
            if i == 0: epsilon = formatted_dos[:,0] * const.AUTOEV          # in eV
            dos.append(formatted_dos[:,1]) 

        if spin_polarized: 
            assert np.mod(ndos,2) == 0, "This is not a spin-polarised calculation"
            ndos = ndos//2
            
        dos = np.asarray(dos).reshape(-1, ndos, nepsilon).transpose(0, 2, 1)
        return epsilon, dos 