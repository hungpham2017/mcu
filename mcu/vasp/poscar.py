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
from ..cell import spg_wrapper, cell_io, cell
from ..cell import utils as cell_utils
from . import utils, vasp_io
import matplotlib as mpl
import matplotlib.pyplot as plt

            
class main(cell.main):
    def __init__(self, poscar='POSCAR'):
        '''Get POSCAR file and return a POSCAR object '''
        self.poscar = vasp_io.POSCAR(poscar)
        self.cell_recip = self.poscar.cell[1]
        self.cell  = utils.cell_to_spgcell(self.poscar.cell, self.poscar.atom)
        self.cell_type = [None, None]
        
    def get_2D_kmesh(self, origin=[0,0,0], krange=[0.1,0.1], plane='xy', npoint=[11,11]):
        '''Get a rectangular k-mesh around a k-point on a plane
        Attribute:
            origin      :   k-point (fractional, reciprocal) coordinate considered at the center of the rectangular
            krange      :   list of the k-mesh window for the 1st and 2sd axis, 0.1 means [-0.1,0.1] in A-1 unit
            plane       :   the plane consider
            npoint      :   number of k-point along the 1st and 2sd axis
        '''
        
        deltax = 2*krange[0]/(npoint[0]-1)
        deltay = 2*krange[1]/(npoint[1]-1)
        X, Y = np.mgrid[-krange[0]:(krange[0]+deltax*0.5):deltax,-krange[1]:(krange[1]+deltay*0.5):deltay]
        coor = np.empty([X.size,3])
        if plane == 'xy':
            coor[:,0] = X.flatten()
            coor[:,1] = Y.flatten()
            coor[:,2] = np.zeros(X.size)
        elif plane == 'xz':
            coor[:,0] = X.flatten()
            coor[:,1] = np.zeros(X.size)   
            coor[:,2] = Y.flatten()     
        elif plane == 'yz':        
            coor[:,0] = np.zeros(X.size) 
            coor[:,1] = X.flatten()
            coor[:,2] = Y.flatten()         

        # To fractional and move to the origin 
        frac_coor = coor.dot(np.linalg.inv(self.cell_recip)) + origin                    
        self.kmesh_2D = frac_coor.dot(self.cell_recip) / (2*np.pi)        # in 2pi/A unit and can be compared with k list in OURCAR
        
        #Write KPOINTS file
        with open('KPOINTS', 'w') as f:
            f.write('Generated mesh by mcu: %s %7.4f %7.4f %2d %2d\n' % (plane, krange[0], krange[1], npoint[0], npoint[1]))	
            f.write('   %4d # + of k-points from IBZKPT\n' % (frac_coor.shape[0]))	    
            f.write('Reciprocal lattice\n')
            f.write('   #This is where you add the k-mesh copied from IBZKPT file, this k-mesh is used to run SCF\n')
            for k in range(frac_coor.shape[0]):
                f.write('%14.10f  %14.10f  %14.10f     %2d\n' % (frac_coor[k,0], frac_coor[k,1], frac_coor[k,2], 0))
                
