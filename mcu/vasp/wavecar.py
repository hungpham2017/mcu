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
import scipy.linalg
from ..utils.misc import check_exist
from . import utils, vasp_io, const
from scipy.fftpack import ifftn
from scipy.io import FortranFile
            
class main:
    def __init__(self, file="WAVECAR", lsorbit=False, vasprun=None):
        '''
           === Attributes ===
            file                : wfs file name, default: WAVECAR
            lsorbit             : collinear or non-colinear calculation
            vasprun             : a vasprun object
        '''
        
        assert check_exist(file), "Cannot find the WAVECAR file. Check the path:" + file
        self._wavecar = open(file, 'r')
        self.lsorbit = lsorbit
        self.ngrid = None
        self.kpts_weight = None
        if vasprun is not None:
            self.lsorbit = vasprun.soc
            self.ispin = vasprun.ispin          # used to check if vasprun.xml and WAVECAR matches
            self.ngrid = vasprun.ngrid 
            self.kpts_weight = vasprun.kpts_weight
            self.get_kpts_nosym = vasprun.get_kpts_nosym
            self.kmesh = vasprun.kmesh
            self.isym = vasprun.isym
            self.isym_symprec = vasprun.isym_symprec
            
        self.read_header()
        self.get_band()

    def read_header(self):
        '''Reading haeder + calc info'''
        self._wavecar.seek(0)
        self._recl, self.nspin, self._rtag = np.array(np.fromfile(self._wavecar, dtype=np.float64, count=3),dtype=int)
        if self._rtag == 45200:
            self.prec = np.complex64
        elif self._rtag == 45210:
            self.prec = np.complex128
        else:
            raise ValueError("Invalid TAG values: {}".format(self._rtag))
            
        # Get kpts, bands, encut, cell info
        self._wavecar.seek(self._recl)
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
        if self.ngrid is None:
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
            for kpt in range(self.nkpts):
                # Read eigenvalues + occ
                rec = 2 + spin*self.nkpts*(self.nbands + 1) + kpt*(self.nbands + 1)
                self._wavecar.seek(rec * self._recl)
                dump = np.fromfile(self._wavecar, dtype=np.float64, count=4+3*self.nbands)
                if spin == 0:
                    self.nplws[kpt] = int(dump[0])
                    self.kpts[kpt] = dump[1:4]
                dump = dump[4:].reshape((-1, 3))
                self.band[spin,kpt,:] = dump[:,0]
                self.co_occ[spin,kpt,:] = dump[:,2]
            
    def get_coeff(self, spin=0, kpt=0, band_list=None):
        '''Extract plw coefficients of the wfn''' 
        
        if band_list is None: band_list = np.arange(self.nbands)
        cg = []
        for band in band_list:
            rec = 3 + spin*self.nkpts*(self.nbands + 1) + kpt*(self.nbands + 1) + band 
            self._wavecar.seek(rec * self._recl)
            dump = np.fromfile(self._wavecar, dtype=self.prec, count=self.nplws[kpt])
            cg.append(dump)
                
        return np.asarray(cg)
        
    def get_gvec(self, kpt=0):
        '''
        Generate the G-vectors that satisfies the following relation
            (G + k)**2 / 2 < ENCUT
        '''

        assert 0 <= kpt < self.nkpts,  'Invalid kpoint index!'
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
        if self.lsorbit: n = 2         # the No. of plw is two times larger for a SOC wfn (up + down) 
        assert Gvec.shape[0] == self.nplws[kpt] / n, 'No. of planewaves not consistent! %d %d %d' % \
                (Gvec.shape[0], self.nplws[kpt], np.prod(self.ngrid))
                         
        return Gvec
        
    def get_unk(self, spin=0, kpt=0, band_list=None, Gp=[0,0,0], ngrid=None, norm=False):
        '''
        Obtain the irreducible pseudo periodic parts of the Bloch function in real space

        Attributes:
            spin    : spin index of the desired KS states, starting from 0
            kpt     : k-point index of the desired KS states, starting from 0
            Gp      : shift the G vectors by Gp 
            gvec    : the G-vectors correspond to the plane-wave coefficients
            Cg      : the plane-wave coefficients. If None, read from WAVECAR
            ngrid   : the FFT grid size
            norm    : thr normalized factor, False: 1/sqrt(N), True: < unk | unk > = 1  
            
            norm=False : the unk will be identical to VASP UNK files
            
             
            u_{n}^{k+Gp}(r) = 1/ norm * \sum_G u_{n}^{k}(G).e^{i(G-Gp)r} 
            
            special case when Gp = 0: 
                u_{n}^{k}(r) = 1/ norm * \sum_G u_{n}^{k}(G).e^{iGr} 
                u_{n}^{k}(r) = FFT(u_{n}^{k}(G)) 
                
        return:
            unk [band-th, nx, ny, nz]
        '''

        if ngrid is None:
            ngrid = np.array(self.ngrid.copy())
        else:
            ngrid = np.array(ngrid, dtype=np.int64)
            assert ngrid.shape[0] == 3, 'Wrong syntax for ngrid'
            assert np.alltrue(ngrid >= self.ngrid), "Minium FT grid size: (%d, %d, %d)" % \
                    (self.ngrid[0], self.ngrid[1], self.ngrid[2])                            

        # The FFT normalization factor
        # the iFFT has a factor 1/N_G, unk exported by VASP does not have this factor
        nband = len(band_list)
        Gp = np.int64(Gp)
        gvec = self.get_gvec(kpt) - Gp
        unk_G = np.zeros([nband, ngrid[0], ngrid[1], ngrid[2]], dtype=self.prec)
        gvec %= ngrid[np.newaxis,:]
        nx, ny, nz = gvec[:,0], gvec[:,1], gvec[:,2]

        if self.lsorbit:
            wfc_spinor = []
            Cg = self.get_coeff(spin, kpt, band_list)
            nplw = Cg.shape[1] // 2
            
            # spinor up
            unk_G[:, nx, ny, nz] = Cg[:, :nplw]
            unk_up = ifftn(unk_G, axes=[1,2,3])
            
            # spinor down
            unk_G[:, nx, ny, nz] = Cg[:, nplw:]
            unk_down = ifftn(unk_G, axes=[1,2,3])

            del Cg, unk_G
            unk = np.hstack([unk_up, unk_down])     # dimension: (nband, nx*2, ny, nz)
            del unk_up, unk_down
        else:
            unk_G[:, nx, ny, nz] = self.get_coeff(spin, kpt)
            unk = ifftn(unk_G, axes=[1,2,3])

        if norm:
            norm = np.einsum('ixyz,jxyz->ij', unk, unk.conj())
            inv_sqrt_norm = scipy.linalg.inv(scipy.linalg.sqrtm(norm))
            norm_unk = np.einsum('ij,jxyz->ixyz', inv_sqrt_norm, unk)
            return norm_unk
        else:
            # Note: ifftn has a norm factor of 1/N, but VASP doesn't have
            return unk * np.prod(ngrid)  
                   
    def get_unk_kpts(self, spin=0, Gp=[0,0,0], ngrid=None, norm=False):
        '''
        Obtain the irreducible pseudo periodic parts of the Bloch function in real space
        at all irreducible kpts
        '''
        unk_kpts = []
        for kpt in range(self.nkpts):
            unk_kpt = self.get_unk(spin=spin, kpt=kpt, Gp=Gp, ngrid=ngrid, norm=norm)
            unk_kpts.append(unk_kpt)
        return unk_kpts    
          
    def get_wave_nosym(self, spin=0, Gp=[0,0,0], ngrid=None, norm=False, match_vasp_kpts=True):
        '''
        This function was meant to collet the entire wave function.
        In practice, it costs huge amount of memory, hence need a lot of work here        
        
        Obtain:
            - the reducible pseudo periodic parts of the Bloch function in real space
            - the reducible band 
        Time-reversal or spatial symmetry is disregarded
        
        TODO:
        match_vasp_kpts is used to debug now with ISYM=-1 and ISYM=0,
        other ISYM, it is not straightforward to match kpts list by VASP
        '''
        
        irred_kpts = self.kpts
        irred_band = self.band[spin]
        irred_unk_kpts = self.get_unk_kpts(spin=spin, Gp=Gp, ngrid=ngrid, norm=norm)

        if self.isym == -1: # No symmetry at all
            self.kpts_nosym = irred_kpts
            self.idx_sym = np.arange(irred_kpts.shape[0])
            return irred_kpts, irred_band, irred_unk_kpts
        elif self.isym == 0: # only time-reversal symmetry
            raise ValueError('Spatial symmetry is not supported yet. ISYM must be -1')
            mapping, kpts_grid = self.get_kpts_nosym(self.isym_symprec, no_spatial=True)
        else: # spatial symmetry
            raise ValueError('Spatial symmetry is not supported yet. ISYM must be -1')
            mapping, kpts_grid = self.get_kpts_nosym(self.isym_symprec)
            
        idx_irred, idx, idx_inverse, idx_count = np.unique(mapping, 
        return_index=True, return_inverse=True, return_counts=True)
 
        kpts = kpts_grid/self.kmesh
        nkpts = kpts.shape[0]
        
        # Check whether the irreducible kpts list matchs the VASP list
        assert np.linalg.norm(kpts[idx_irred] - irred_kpts) < 1e-6, \
                    "The calculated irreducible kpts list does not match the one used in VASP"
        assert np.linalg.norm(idx_count/nkpts - self.kpts_weight[:,0]) < 1e-6, \
                    "The calculated kpts weight list does not match the one used in VASP"
        assert abs(self.nspin - self.ispin) < 1e-8, \
                    "The ispin in WAVECAR and vasprun.xml are not the same. Check if they are from the same calculation"
        
        if match_vasp_kpts:
            nirred = len(idx_irred)
            orig_item = []
            image_item = []
            temp = np.arange(kpts.shape[0])
            for i in range(nirred):
                irred_item = temp[idx_inverse==i].tolist()
                orig_item.append(irred_item[0])
                image_item.append(irred_item[1:])
            new_idx = orig_item + sum(image_item,[])
            kpts = kpts[new_idx]
            idx_inverse = idx_inverse[new_idx]
            
        band = []
        unk = []
        for kpt in range(nkpts):
            band.append(irred_band[idx_inverse[kpt]])
            k_reserv = abs(kpts[kpt] + self.kpts[idx_inverse[kpt]]).sum() < 1e-7
            if abs(kpts[kpt]).sum() > 1e-7 and k_reserv:
                unk_kpt = irred_unk_kpts[idx_inverse[kpt]].conj()
            else:
                unk_kpt = irred_unk_kpts[idx_inverse[kpt]]
            unk.append(unk_kpt)

        self.kpts_nosym = kpts
        self.idx_sym = idx_inverse
        
        return kpts, np.asarray(band), unk
        
    def get_unk_kpt(self, spin=0, kpt=0, band_list=None, Gp=[0,0,0], ngrid=None, norm=False): 
        '''
        Obtain the irreducible pseudo periodic parts of the Bloch function in real space
        at all abitrary kpts considering symmetry
        Must run after calling get_wave_nosym
        '''

        assert 0 <= kpt < self.kpts.shape[0], "Invalid value of kpt"
        unk_kpt = self.get_unk(spin=spin, kpt=kpt, band_list=band_list, Gp=Gp, ngrid=ngrid, norm=norm)
        return unk_kpt

        
    def write_vesta(self, unk, realonly=False, poscar='POSCAR', filename='unk', ncol=10):
        '''
        Save the real/imag space pseudo-wavefunction as vesta format.
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
            
    def export_unk(self, spin=0, ngrid=None, only_irred=False):
        '''
        Export the periodic part of BF in a real space grid for plotting with wannier90
        '''    
        
        if self.lsorbit:
            spin_str = '.NC'
            spin = 0
        elif spin == 0:
            spin_str = '.1'
        else:
            spin_str = '.2'
            
        if ngrid is None:
            ngrid = self.ngrid.copy()
        else:
            ngrid = np.array(ngrid, dtype=np.int64)
            assert ngrid.shape[0] == 3, 'Wrong syntax for ngrid'
            assert np.alltrue(ngrid >= self.ngrid), "Minium FT grid size: (%d, %d, %d)" % \
                    (self.ngrid[0], self.ngrid[1], self.ngrid[2])

        if only_irred:
            nkpts = self.nkpts
            unk_list = self.get_unk_kpts(spin, ngrid=ngrid)
        else:
            kpts, band, unk_list = self.get_wave_nosym(spin, ngrid=ngrid)
            nkpts = kpts.shape[0]
        
        if self.lsorbit:
            for kpt in range(nkpts):
                unk_file = FortranFile('UNK' + "%05d" % (kpt + 1) + spin_str, 'w')
                unk_file.write_record(np.asarray([ngrid[0], ngrid[1], ngrid[2], kpt + 1, self.nbands], dtype=np.int32)) 
                for band in range(self.nbands):
                    unk_up = unk_list[kpt][band,:ngrid[0],:,:].T.flatten()
                    unk_down = unk_list[kpt][band,ngrid[0]:,:,:].T.flatten()
                    unk_file.write_record(unk_up)
                    unk_file.write_record(unk_down)                    
                unk_file.close()
        else:
            for kpt in range(nkpts):
                unk_file = FortranFile('UNK' + "%05d" % (kpt + 1) + spin_str, 'w')
                unk_file.write_record(np.asarray([ngrid[0], ngrid[1], ngrid[2], kpt + 1, self.nbands], dtype=np.int32)) 
                for band in range(self.nbands):
                    unk = unk_list[kpt][band].T.flatten()
                    unk_file.write_record(unk)                    
                unk_file.close()