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
           
############ CRYSTAL output ##################
primitive_cell_MATCH = re.compile(r'''
[\w\W]*
 LATTICE [ ]* PARAMETERS [ ]* \( [ ]* ANGSTROMS [ ]* AND [ ]* DEGREES [ ]* \) [ ]* - [ ]* BOHR [ ]* = [ ]* (?P<BOHR>\S+) [ ]* ANGSTROM
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n\n)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

lattice_MATCH = re.compile(r'''
[\w\W]*
  A [ ]* B [ ]* C [ ]* ALPHA [ ]* BETA [ ]* GAMMA
(?P<content>
  [\s\S]*?(?=\n.*?[ ] \*|$)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

irred_atoms_MATCH = re.compile(r'''
[\w\W]*
 ATOM [ ]* (?P<X>\S+) [ ]* (?P<Y>\S+) [ ]* (?P<Z>\S+)
 [\*]*
 (?P<content>
  [\w\W]*
)
''', re.VERBOSE)

atoms_MATCH = re.compile(r'''
[ ]*
 (?P<order>\d+) [ ]* (?P<asym>\S+) [ ]* (?P<Z>\d+) [ ]* (?P<atom>\S+) [ ]* (?P<x>\S+) [ ]* (?P<y>\S+) [ ]*  (?P<z>\S+)
''', re.VERBOSE)

nshells_MATCH = re.compile(r'''
[\w\W]* NUMBER [ ]* OF [ ]* SHELLS [ ]* (?P<nshells>\d+)
''', re.VERBOSE)

nao_MATCH = re.compile(r'''
[\w\W]* NUMBER [ ]* OF [ ]* AO [ ]* (?P<nao>\d+)
''', re.VERBOSE)

nelec_MATCH = re.compile(r'''
[\w\W]*  N. [ ]* OF [ ]* ELECTRONS [ ]* PER [ ]* CELL [ ]* (?P<nelec>\d+)
''', re.VERBOSE)

nelec_core_MATCH = re.compile(r'''
[\w\W]* CORE [ ]* ELECTRONS [ ]* PER [ ]* CELL [ ]* (?P<nelec>\d+)
''', re.VERBOSE)
 
kpts_MATCH = re.compile(r'''
[\w\W]*
 [\*]* [ ]* K [ ]* POINTS [ ]* COORDINATES [ ]* \( [ ]* OBLIQUE [ ]* COORDINATES [ ]* IN [ ]* UNITS [ ]* OF [ ]* IS [ ]* = [ ]* (?P<IS>\d+) [ ]* \)
 (?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n\n)
)
''', re.VERBOSE)

kpt_MATCH = re.compile(r'''
[ ]*
 (?P<order>\d+)-(?P<label>\S+)\( [ ]* (?P<kx>\S+) [ ]* (?P<ky>\S+) [ ]* (?P<kz>\S+) [ ]*\)
''', re.VERBOSE)

basisset_MATCH = re.compile(r'''
[\w\W]*
   ATOM [ ]* X\(AU\) [ ]* Y\(AU\) [ ]* Z\(AU\) [ ]* N. [ ]* TYPE [ ]* EXPONENT [ ]* S [ ]* COEF [ ]* P [ ]* COEF [ ]* D/F/G [ ]* COEF
 [\*]*
 (?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n\n)
)
''', re.VERBOSE)

basis_atom_MATCH = re.compile(r'''
[ ]+
 (?P<order>\d+) [ ]+ (?P<atom>\S+) [ ]+ (?P<x>\S+) [ ]+ (?P<y>\S+) [ ]+ (?P<z>\S+)XX
 (?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n))
''', re.VERBOSE)
   
basis_function_MATCH = re.compile(r'''
[ ]+
 (?P<order>\d+) [ ]+ (?P<orb>\S+) [ ]+ XX
''', re.VERBOSE)
   

def read_out(filename):
    '''Read prefix.out file'''
 
    with open(filename, 'r') as data_file:
        data = data_file.read()
        
        # Get the primitive lattice
        primitive_cell = primitive_cell_MATCH.match(data).group('content')
        temp = lattice_MATCH.match(primitive_cell).group('content')
        lattice = np.float64(temp.split())
        lattice = cell_utils.convert_lattice(lattice)
                
        # Get irreducible atoms
        irred_atoms = irred_atoms_MATCH.match(primitive_cell).group('content')        
        atom = []
        atom_position = []
        for atom_data in atoms_MATCH.finditer(irred_atoms):
            content = atom_data.groupdict()
            atom.append(content['atom'].capitalize())   
            atom_position.append([content['x'], content['y'], content['z']])  

        # electronic information
        nshells = int(nshells_MATCH.match(data)['nshells'])
        nao = int(nao_MATCH.match(data)['nao'])
        nelec = int(nelec_MATCH.match(data)['nelec'])
        nelec_core = int(nelec_core_MATCH.match(data)['nelec'])

        # Get kpts:
        temp = kpts_MATCH.match(data)
        IS = int(temp['IS'])
        kpts_data = temp.group('content')
        kpts_label = []
        kpts = []
        for kpt_data in kpt_MATCH.finditer(kpts_data):
            content = kpt_data.groupdict()
            kpts_label.append(content['label']) 
            kpts.append([content['kx'], content['ky'], content['kz']])  
        kpts = np.float64(kpts) / IS
        
        # Get basis function
        basisset = basisset_MATCH.match(data).group('content')
        basis_function = basisset.replace('\n                                         ','XX')
        basis_atom = basisset.replace('\n              ','XX')

        basis = {}
        for irred_atom in basis_atom_MATCH.finditer(basis_atom):
            content = irred_atom.groupdict()
            element = content['atom'].capitalize()
            temp = content['content'].replace('\n                                         ','XX')
            function = []
            for orb in basis_function_MATCH.finditer(temp):
                function.append(orb.groupdict()['orb'].lower())
            basis[element] = function
        
        out = {}
        out['nelec'] = nelec
        out['nelec core'] = nelec_core
        out['nao'] = nao
        out['basis'] = basis 
        out['nshells'] = nshells      
        out['lattice'] = lattice
        out['atom'] = atom
        out['atom_position'] = np.float64(atom_position)
        out['kpts'] = kpts
        
        
    return out
    
        
def read_BAND(filename):
    '''Read prefix.BAND/BAND.DAT file'''
    
    if not check_exist(filename):
        if check_exist("BAND.DAT"):
            print('Cannot find : ' + filename, '. BAND.DAT is used instead')
            filename = "BAND.DAT"
        else:
            assert False, 'Cannot find the BAND file. Check the path: ' + filename
    
    with open(filename, 'r') as data_file:
        data = data_file.readlines()
        nkpt, nband, nspin = np.int64([data[0].split()[i] for i in [2,4,6]])
        npath = int(data[1].split()[2])
        start = npath + 3
        block_length = 19 + 2*npath + nkpt
        band = []
        for spin in range(nspin):
            block = data[(start + block_length*spin):(start + block_length*(spin + 1))]
            symm_k_coor = np.float64([block[14 + 2*i].split()[-1] for i in range(npath + 1)])
            band_block = block[(18 + 2*npath):-1]
            band_spin = np.float64("\n".join(band_block).split())
            band_spin = band_spin.reshape(nkpt, nband+1)
            band.append(band_spin[:,1:])
            kpath_abs = band_spin[:,0] 
            efermi = np.float64(block[-1].split()[-1])

        band = np.float64(band)
            
    return symm_k_coor, kpath_abs, band, efermi

def read_DOSS(filename):
    '''Read prefix.DOSS/DOSS.DAT file'''
    
    if not check_exist(filename):
        if check_exist("DOSS.DAT"):
            print('Cannot find : ' + filename, '. DOSS.DAT is used instead')
            filename = "DOSS.DAT"
        else:
            assert False, 'Cannot find the DOSS file. Check the path: ' + filename

    with open(filename, 'r') as data_file:
        data = data_file.readlines()
        nepts, nproj, nspin = np.int64([data[0].split()[i] for i in [2,4,6]]) 
        block_length = nepts + 2
        
        dos = []
        for spin in range(nspin):
            block = data[(3 + block_length*spin):(3 + block_length*(spin + 1))]
            dos_block = block[1:-1]
            dos_spin = np.float64("".join(dos_block).split())
            dos_spin = dos_spin.reshape(nepts, nproj + 1)
            dos.append(dos_spin)
            epts = dos_spin[:,0] 
            efermi = np.float64(block[-1].split()[-1])
            
        dos = np.float64(dos)
    
    return epts, dos, efermi
        

f25_BAND_MATCH = re.compile(r'''
[ ]*
 -\%-(?P<ihferm>\d)BAND [ ]* (?P<nband>\d+) [ ]* (?P<nkp>\d+)(?P<dum>[\s\S]{12}?(?=.[ ]|))(?P<dk>[\s\S]{12}?(?=.[ ]|))(?P<efermi>[\s\S]{12}?(?=.[ ]|))
 \n(?P<emin>[\s\S]{12}?(?=.[ ]|))(?P<emax>[\s\S]{12}?(?=.[ ]|))\n
 [ ]* (?P<val1>\d+) [ ]* (?P<val2>\d+) [ ]* (?P<val3>\d+) [ ]* (?P<val4>\d+) [ ]* (?P<val5>\d+) [ ]* (?P<val6>\d+) 
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n-%-)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

f25_DOSS_MATCH = re.compile(r'''
[ ]*
 -\%-(?P<ihferm>\d)DOSS [ ]* (?P<nrow>\d+) [ ]* (?P<ncol>\d+)(?P<dx>[\s\S]{12}?(?=.[ ]|))(?P<dy>[\s\S]{12}?(?=.[ ]|))(?P<efermi>[\s\S]{12}?(?=.[ ]|))
 \n(?P<emin>[\s\S]{12}?(?=.[ ]|))(?P<emax>[\s\S]{12}?(?=.[ ]|))\n
 [ ]* (?P<val1>\d+) [ ]* (?P<val2>\d+) [ ]* (?P<val3>\d+) [ ]* (?P<val4>\d+) [ ]* (?P<val5>\d+) [ ]* (?P<val6>\d+) 
(?P<content>
  [\s\S]*?(?=\n.*?[ ] $|\n-%-)  # match everything until next blank line or EOL
)
''', re.VERBOSE)

        
def read_f25(filename):
    '''Read prefix.f25/fort.25 file'''
    with open(filename, 'r') as data_file:
        data = data_file.read()
        
        # Get BAND
        if f25_BAND_MATCH.match(data + '\n-%-') is None:
            out_BAND = None
        else:
            out_BAND = []
            for block in f25_BAND_MATCH.finditer(data + '\n-%-'):
                ihferm = int(block['ihferm'])
                nband = int(block['nband'])
                nkp = int(block['nkp'])
                dum = float(block['dum'])
                dk = float(block['dk'])
                efermi = float(block['efermi'])
                emin = float(block['emin'])
                emax = float(block['emax'])
                values = np.int64([block['val1'], block['val2'], block['val3'], block['val4'], block['val5'], block['val6']])
                band_data = block.group('content').replace("\n", "")
                eigenvals = np.float64(list(map(''.join, zip(*[iter(band_data)]*12))))
                
                block_data = {}
                block_data['ihferm'] = ihferm
                block_data['nband'] = nband
                block_data['nkp'] = nkp
                block_data['dum'] = dum
                block_data['dk'] = dk  
                block_data['efermi'] = efermi
                block_data['emin'] = emin
                block_data['emax'] = emax
                block_data['values'] = values    
                block_data['eigenvals'] = eigenvals
                out_BAND.append(block_data)
                
        # Get DOSS
        if f25_DOSS_MATCH.match(data + '\n-%-') is None:
            out_DOSS = None
        else:
            out_DOSS = []
            for block in f25_DOSS_MATCH.finditer(data + '\n-%-'):
                ihferm = int(block['ihferm'])
                nrow = int(block['nrow'])
                ncol = int(block['ncol'])
                dx = float(block['dx'])
                dy = float(block['dy'])
                efermi = float(block['efermi'])
                emin = float(block['emin'])
                emax = float(block['emax'])
                values = np.int64([block['val1'], block['val2'], block['val3'], block['val4'], block['val5'], block['val6']])
                band_data = block.group('content').replace("\n", "")
                dos = np.float64(list(map(''.join, zip(*[iter(band_data)]*12))))
                
                block_data = {}
                block_data['ihferm'] = ihferm
                block_data['nrow'] = nrow
                block_data['ncol'] = ncol
                block_data['dx'] = dx
                block_data['dy'] = dy 
                block_data['efermi'] = efermi
                block_data['emin'] = emin
                block_data['emax'] = emax
                block_data['values'] = values    
                block_data['dos'] = dos
                out_DOSS.append(block_data)

    return out_BAND, out_DOSS
       