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

'''This module is modied from the vaspwfc.py from QijingZheng's project.
    ref: https://github.com/QijingZheng/VaspBandUnfolding/blob/master/vaspwfc.py)
'''

import numpy as np
from mcu.vasp import utils, vasp_io, const
from scipy.fftpack import fftfreq, fftn, ifftn
            
class main:
    def __init__(self, file="WAVECAR", lsorbit=False):
        '''The wavecar manipulation is modied from the vaspwfc.py from QijingZheng's project
            ref: https://github.com/QijingZheng/VaspBandUnfolding/blob/master/vaspwfc.py)
        '''
        if not utils.check_exist(file):
            print('Cannot find the WAVECAR file. Check the path:', file)
            self.success = False
        else: 
            self._wavecar = open(file, 'rb')
            self._lsorbit = lsorbit
            self.success = True
            self.read_header()
            self.get_band()

    def read_header(self):
        '''Reading haeder + calc info'''
        self._wavecar.seek(0)
        self.recl, self.nspin, self.rtag = np.array(np.fromfile(self._wavecar, dtype=np.float64, count=3),dtype=int)
        if self.rtag == 45200:
            self.prec = np.complex64
        elif self.rtag == 45210:
            self.prec = np.complex128
        else:
            raise ValueError("Invalid TAG values: {}".format(self.rtag))
            
        # Get kpts, bands, encut, cell info
        self._wavecar.seek(self.recl)
        dump = np.fromfile(self._wavecar, dtype=np.float64, count=12)
        self.nkpts  = int(dump[0])                              # No. of k-points
        self.nbands = int(dump[1])                              # No. of bands
        self.encut  = dump[2]                                   # Energy cutoff
        lattice  = dump[3:].reshape((3,3))                      # real space supercell basis
        volume  = np.linalg.det(lattice)                        # real space supercell volume, unit: A^-3
        recip_lattice  = 2*np.pi*np.linalg.inv(lattice).T       # reciprocal space supercell volume, unit: A^-1
        self.cell = (lattice, recip_lattice, None, volume)
        
        # Estimating FFT grid size, 
        # old version: 2 * CUTOFF + 1, relax this condition works for the LSORBIT 
        norm = np.linalg.norm(lattice, axis=1)
        CUTOFF = np.ceil(np.sqrt(self.encut/const.RYTOEV) / (2*np.pi/(norm / const.AUTOA)))
        self.ngrid = np.array(2 * CUTOFF + 3, dtype=np.int64)
        
    def get_band(self):
        '''Extract band, occ'''
        
        self.nplws = np.zeros(self.nkpts, dtype=np.int64)
        self.kpts = np.zeros((self.nkpts, 3), dtype=np.float64)
        self.band = np.zeros((self.nspin, self.nkpts, self.nbands), dtype=np.float64)
        self.co_occ  = np.zeros((self.nspin, self.nkpts, self.nbands), dtype=np.float64)           
        for spin in range(self.nspin):
            cg_spin = []
            for kpt in range(self.nkpts):
                # Read eigenvalues + occ
                rec = 2 + spin*self.nkpts*(self.nbands + 1) + kpt*(self.nbands + 1)
                self._wavecar.seek(rec * self.recl)
                dump = np.fromfile(self._wavecar, dtype=np.float64, count=4+3*self.nbands)
                if spin == 0:
                    self.nplws[kpt] = int(dump[0])
                    self.kpts[kpt] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self.band[spin,kpt,:] = dump[:,0]
                self.co_occ[spin,kpt,:] = dump[:,2]
            
    def get_coeff(self, spin=0, kpt=0, norm=False):
        '''Extract plw coefficients of the wfn''' 

        if kpt >= self.nkpts:
            raise ValueError("kpt must be smaller than the maximum index for kpt", self.nkpts)
            
        #TODO: check spin value

        cg = []
        for band in range(self.nbands):
            rec = 3 + spin*self.nkpts*(self.nbands + 1) + kpt*(self.nbands + 1) + band 
            self._wavecar.seek(rec * self.recl)
            dump = np.fromfile(self._wavecar, dtype=self.prec, count=self.nplws[kpt])
            cg.append(np.asarray(dump, dtype=np.complex128))
            
        if norm: cg = cg/np.linalg.norm(cg)

        return np.asarray(cg)
        
    def get_gvec(self, kpt=0):
        '''
        Generate the G-vectors that satisfies the following relation
            (G + k)**2 / 2 < ENCUT
        '''
        assert 0 <= kpt  <= self.nkpts - 1,  'Invalid kpoint index!'
        kvec = self.kpts[kpt]
        
        # Fast algorithm without for loop
        fx = np.hstack([np.arange(self.ngrid[0]//2 + 2), -np.flip(np.arange(1,self.ngrid[0]//2))])
        fy = np.hstack([np.arange(self.ngrid[1]//2 + 2), -np.flip(np.arange(1,self.ngrid[1]//2))])
        fz = np.hstack([np.arange(self.ngrid[2]//2 + 2), -np.flip(np.arange(1,self.ngrid[2]//2))])
        y, z, x = np.meshgrid(fy, fz, fx, indexing='xy')
        kgrid = np.asarray(np.hstack([x.reshape(-1,1),y.reshape(-1,1),z.reshape(-1,1)]))

        # Kinetic_Energy = (G + k)**2 / 2
        # HSQDTM    =  hbar**2/(2*ELECTRON MASS)
        recip_lattice = self.cell[1]
        KENERGY = const.HSQDTM * np.linalg.norm(np.dot(kgrid + kvec[np.newaxis,:] , recip_lattice), axis=1)**2
        
        # find Gvectors where (G + k)**2 / 2 < ENCUT
        Gvec = kgrid[np.where(KENERGY < self.encut)[0]]
        
        # Check if the Gvec is consistent with Gvec generated by VASP
        n = 1
        if self._lsorbit: n = 2         # the No. of plw is two times larger for a SOC wfn (up + down) 
        assert Gvec.shape[0] == self.nplws[kpt] / n, 'No. of planewaves not consistent! %d %d %d' % \
                (Gvec.shape[0], self.nplws[kpt], np.prod(self.ngrid))
                         
        return Gvec
        
    def get_unk(self, spin=0, kpt=1, band=1, ngrid=None, norm_u=True, norm_c=False):
        '''
        Obtain the pseudo periodic part of the Bloch function in real space

        Attributes:
            spin    : spin index of the desired KS states, starting from 1
            kpt     : k-point index of the desired KS states, starting from 1
            band    : band index of the desired KS states, starting from 1
            gvec    : the G-vectors correspond to the plane-wave coefficients
            Cg      : the plane-wave coefficients. If None, read from WAVECAR
            ngrid   : the FFT grid size
            norm_u  : whether to normalzie the u
            norm_c  : whether to normalzie cg            

        The wavefunctions can be normalized such that:

                        \sum_{ijk} | \phi_{ijk} | ^ 2 = 1
            
        '''
        kpt -= 1
        band -= 1

        if ngrid is None:
            ngrid = 2 * self.ngrid.copy()
        else:
            ngrid = np.array(ngrid, dtype=np.int64)
            assert ngrid.shape[0] == 3, 'Wrong syntax for ngrid'
            assert np.alltrue(ngrid >= self.ngrid), "Minium FT grid size: (%d, %d, %d)" % \
                    (self.ngrid[0], self.ngrid[1], self.ngrid[2])

        # default normalization factor so that 
        # \sum_{ijk} | \phi_{ijk} | ^ 2 = 1
        normfac = np.sqrt(np.prod(ngrid)) if norm_u else 1.0

        gvec = self.get_gvec(kpt)
        unk = np.zeros(ngrid, dtype=np.complex128)
        gvec %= ngrid[np.newaxis,:]
        nx, ny, nz = gvec[:,0], gvec[:,1], gvec[:,2]

        if self._lsorbit:
            wfc_spinor = []
            Cg = self.get_coeff(spin, kpt, norm_c)[band]  
            nplw = Cg.shape[0] // 2
            
            # spinor up
            unk[nx, ny, nz] = Cg[:nplw]
            wfc_spinor.append(ifftn(unk))
            
            # spinor down
            unk[:,:,:] = 0.0j
            unk[nx, ny, nz] = Cg[nplw:]
            wfc_spinor.append(ifftn(unk))

            del Cg
            return np.asarray(wfc_spinor)*normfac
            
        else:
            unk[nx, ny, nz] = self.get_coeff(spin, kpt, norm_c)[band] 
            return ifftn(unk * normfac)
                
    def write_vesta(self, unk, realonly=False, poscar='POSCAR', filename='unk',
                   ncol=10):
        '''
        Save the real space pseudo-wavefunction as vesta format.
        '''
        nx, ny, nz = unk.shape
        try:
            pos = open(poscar, 'r')
            head = ''
            for line in pos:
                if line.strip():
                    head += line
                else:
                    break
            head += '\n%5d%5d%5d\n' % (nx, ny, nz)
        except:
            raise IOError('Failed to open %s' % poscar)

        # Faster IO
        nrow = unk.size // ncol
        nrem = unk.size % ncol
        fmt = "%16.8E"

        psi = unk.copy()
        psi = psi.flatten(order='F')
        psi_h = psi[:nrow * ncol].reshape((nrow, ncol))
        psi_r = psi[nrow * ncol:]

        # Write the real part
        with open(filename + '_r.vasp', 'w') as out:
            out.write(head)
            out.write(
                '\n'.join([''.join([fmt % xx for xx in row])
                           for row in psi_h.real])
            )
            out.write("\n" + ''.join([fmt % xx for xx in psi_r.real]))
            
        # Write the imaginary part            
        if not realonly:
            with open(filename + '_i.vasp', 'w') as out:
                out.write(head)
                out.write(
                    '\n'.join([''.join([fmt % xx for xx in row])
                               for row in psi_h.imag])
                )
                out.write("\n" + ''.join([fmt % xx for xx in psi_r.imag]))    
                