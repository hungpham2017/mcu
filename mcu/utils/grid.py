#!/usr/bin/env python -u 
'''
pDMET: Density Matrix Embedding theory for Periodic Systems
Copyright (C) 2018 Hung Q. Pham. All Rights Reserved.
A few functions in pDMET are modifed from QC-DMET Copyright (C) 2015 Sebastian Wouters

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

import os, sys, datetime
import numpy as np        


def cartesian_prod(arrays, out=None, order='C'):
    '''
    Generate lattice vectors
    '''
    arrays = [np.asarray(x) for x in arrays]
    dtype = np.result_type(*arrays)
    nd = len(arrays)
    dims = [nd] + [len(x) for x in arrays]
    if out is None:
        out = np.empty(dims, dtype)
    else:
        out = np.ndarray(dims, dtype, buffer=out)
    tout = out.reshape(dims)
    shape = [-1] + [1] * nd
    for i, arr in enumerate(arrays):
        tout[i] = arr.reshape(shape[:nd-i])
    return tout.reshape((nd,-1),order=order).T
    
def periodic_grid(lattice, grid=[30,30,30], supercell=[1,1,1], order='C'):
	'''
	Generate a periodic grid for the unit/computational cell in F/C order
    Note: coords has the same unit as lattice
	'''	
	ngrid = np.asarray(grid)
	qv = cartesian_prod([np.arange(-ngrid[i]*(supercell[i]//2),ngrid[i]*((supercell[i]+1)//2)) for i in range(3)], order=order)   
	a_frac = np.einsum('i,ij->ij', 1./ngrid, lattice)
	coords = np.dot(qv, a_frac)
    
	# Compute weight    
	ngrids = np.prod(grid)
	ncells = np.prod(supercell)
	weights = np.empty(ngrids*ncells)
	vol = abs(np.linalg.det(lattice))
	weights[:] = vol / ngrids / ncells   
	return coords, weights