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

import os
import numpy as np
from . import wavecar
from ..utils.grid import periodic_grid


def get_parity(vasprun, spin=0, kpt=0, band_list=None, time_reversal=True):
    '''Get the parity of the Bloch function for the centrosymmetric crystal
       The inversion center is assumed to be at the center of the unit cell
       == Attributes ==
        kpt             : k-point, start at 0
        band_list       : list of bands to be computed, start from 1 to match with VASP's convention
                          by default, all the occupied band will be computed
        time_reversal   : if True and the wave function is non-collinear, alpha and beta band will be degenerate
                           hence, only alpha band will be computed.
                          
       == Return ==
        parity = <psi_{n,k}| P^{hat} |<psi_{n,k}> = <psi_{n,k}|<psi_{n,-k}>
        '''  
    assert kpt < vasprun.kpts.shape[0], "kpt must be less than " + str(vasprun.kpts.shape[0])
    vasprun.isym = -1
    wave = wavecar.main(vasprun=vasprun)
    nelec = vasprun.nelec
    nbands = vasprun.nbands
    soc = vasprun.soc 
    if soc:
        nocc = int(nelec)
    else:
        nocc = int(nelec/2)
        
    if band_list is None:
        if soc and time_reversal:
            band_list =  2*np.arange(0,nocc//2)
        elif soc and not time_reversal:
            band_list =  np.arange(0,nocc)
        else:
            band_list =  np.arange(0,nocc)
    else:
        assert (np.asarray(band_list) < nbands).all(), "There is only " + str(nocc) + " occupied bands. Check your band list"
    
    lattice = vasprun.cell[0]
    coords, weights = periodic_grid(lattice, wave.ngrid, order = 'F')
    recip_lattice = 2*np.pi*np.linalg.inv(lattice).T
    kpts_abs = vasprun.kpts.dot(recip_lattice)

    delta = 1
    not_symmetric = False
    print('Band    Xi       Parity', flush=True)
    for n in band_list:
        u = wave.get_uk(spin=0, kpt=kpt, band=n,  minus_r=False, norm=True)
        exp = np.exp(1j*coords.dot(kpts_abs[kpt])).reshape(wave.ngrid, order = 'F') 
        if soc: exp = np.vstack([exp,exp])
        psi = np.einsum('xyz,xyz->xyz', exp, u)
        u_minus = wave.get_uk(spin=0, kpt=kpt, band=n,  minus_r=True, norm=True)
        exp = np.exp(1j*(np.ones(3).dot(lattice) - coords).dot(kpts_abs[kpt])).reshape(wave.ngrid, order = 'F') 
        if soc: exp = np.vstack([exp,exp])
        psi_minus = np.einsum('xyz,xyz->xyz', exp, u_minus)
        P = (psi.conj() * psi_minus).sum()
        if abs(P) < 0.9: not_symmetric = True
        sign = int(np.sign(P.real))
        sign_ = "+"
        if sign==-1:
            sign_ = "-"
        print('{0:3d}     {1:3s}    {2:3.4f}'.format(n + 1, sign_, P.real), flush=True)
        delta *= sign
    print("delta = ", delta, flush=True)
    if not_symmetric: print("WARNING: the Bloch function is not symmetric around the inversion center! The eigenvalues of the parity operator is not 1.")