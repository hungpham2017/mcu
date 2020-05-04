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
Original from CRYSTAL17 manual
0 S 1 s 
1 SP 4 s, x, y, z 
2 P 3 x, y, z 
3 D 5 2z2 − x2 − y2, xz, yz, x2 − y2, xy 
4 F 7 (2z2 − 3x2 − 3y2)z, (4z2 −x2 − y2)x, (4z2 − x2 − y2)y, (x2 − y2)z, xyz, (x2 − 3y2)x, (3x2 − y2)y
'''

'''
Modied to fit 'standard' (simpler) notation
S     1 s 
P 3 x,y,z 
D 5 dz2, dxz, dyz, dx2−y2, dxy 
F 7 dz3, fxz2, fyz2, fz(x2−y2), fxyz, fx(x2−3y2), fy(3x2−y2)
'''

lm_basis_dict = {'s':['s'], 'p':['px','py','pz'], 'd':['dz2', 'dxz','dyz','dx2-y2','dxy'],
                 'f':['fz3', 'fxz2', 'fyz2', 'fz(x2-y2)', 'fxyz', 'fx(x2-3y2)', 'fy(3x2-y2)']}
lm_basis_list = sum(list(lm_basis_dict.values()),[])

def basis_short2long(basis):
    '''Give a short list of basis set, return the full list'''
    long_basisset = {}
    for element in basis:
        orbs = [lm_basis_dict[ao] for ao in basis[element]]
        orbs = sum(orbs,[])
        long_basisset[element] = orbs
    return long_basisset
    
def orb_index(atom, basis, atom_th, lm):
    '''Give an atom list, a basis set, a lm string return a list of orbital order (start from 1)'''
    
    assert lm in lm_basis_list or lm is None, "The orbital " + lm + " cannot be found in the data base :" + ', '.join(lm_basis_list)
    atom_1storb_idx = []
    count = 0
    for i, atm in enumerate(atom):
        atom_1storb_idx.append(count)
        count += sum([len(lm_basis_dict[orb]) for orb in basis[atm]])
     
    atom_name = atom[atom_th]
    orb_idx = []
    start = atom_1storb_idx[atom_th] + 1
    long_basis = basis_short2long(basis)
    if lm is None:
        orb_idx = [i + start for i, orb in enumerate(long_basis[atom_name])]
    else:
        orb_idx = [i + start for i, orb in enumerate(long_basis[atom_name]) if orb == lm]
    
    return orb_idx