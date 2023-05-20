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
from ..utils import plot, str_format
from ..utils.misc import check_exist
from ..cell import utils as cell_utils
from ..cell import cell
from ..vasp import const
from . import utils as pscf_utils
from . import pscf_io
        
class main(cell.main, plot.main):
    def __init__(self,  pyscf_cell):
        self.pyscf_cell = pyscf_cell
        self.get_info()
        
############ General #################
        
    def get_info(self):    
        '''Extract basis information from the vasprun.xml'''
                
        self.nelec = self.pyscf_cell.nelectron
        self.nao = self.pyscf_cell.nao
        self.nbands = self.nao
        self.soc = False
        
        # Make a cell object in the spglib format
        atom = []
        abs_positions = []
        for atm in self.pyscf_cell.atom:
            atom.append(atm[0])
            abs_positions.append(atm[1])
            
        self.atom  = atom
        self.natom  = len(self.atom)
        self.element = list(dict.fromkeys(self.atom))
        lattice = self.pyscf_cell.lattice_vectors()*const.AUTOA
        positions = abs_positions @ np.linalg.inv(lattice)
        numbers = cell_utils.convert_atomtype(self.atom)
        self.cell = (lattice, positions, numbers)
        self.band = None
        self.kpts = None        
        
############ Plotting ################# 
    def make_kpts(self, kpath, npoint=20):
        '''Giving a 1D k-path, return the kpts list in frac or absolute coordinates'''
        labels, coords = str_format.format_klabel(kpath)
        frac_kpts = []
        for path_th in range(len(coords) - 1):
            temp = np.arange(npoint)/(npoint - 1)
            path = np.ones([npoint,3]) * np.asarray([temp]*3).T * (coords[path_th + 1] - coords[path_th]) + coords[path_th]    
            if path_th != len(coords) - 2:
                path = path[:-1,:]
            frac_kpts.append(path)

        frac_kpts = np.concatenate(frac_kpts)
        lattice = self.cell[0]
        recip_lattice = 2 * np.pi * np.linalg.inv(lattice).T # Get the reciprocal lattice in the row vector format
        abs_kpts = frac_kpts.dot(recip_lattice) * const.AUTOA       # in AU for PySCF 
        
        return frac_kpts, abs_kpts      

    def set_kpts_bands(self, list_or_tuple_or_filename):
        '''Set kpts (fractional) and bands from PySCF calculation'''
        if isinstance(list_or_tuple_or_filename, list) or isinstance(list_or_tuple_or_filename, tuple):
            self.kpts, band = list_or_tuple_or_filename
        elif isinstance(list_or_tuple_or_filename, str):
            self.kpts, band = pscf_io.load_kpts_bands(list_or_tuple_or_filename)
        
        # The number of MO at each kpt may be different, causing an insistent number of bands
        # to avoid inconvenience, only 
        nkpts = self.kpts.shape[0]
        
        nmo_smallest = np.min([len(mo_energy) for mo_energy in band[0]])
        mo_energy_kpts = []
        mo_coeff_kpts = []
        for kpt in range(nkpts):
            mo_energy_kpts.append(band[0][kpt][:nmo_smallest])
            mo_coeff_kpts.append(band[1][kpt][:,:nmo_smallest])
        self.band = [mo_energy_kpts, mo_coeff_kpts]
            
    def save_kpts_bands(self, filename, list_or_tuple_or_filename):
        pscf_io.save_kpts_bands(filename, list_or_tuple_or_filename)
    
    def get_efermi(self, band=None):
        '''E_fermi is assumed to be the valence band maximum
        '''
        if band is None: 
            assert self.band is not None, "You need to provide bands calculated by PySCF kks.get_bands function"
            band = np.asarray(self.band[0])
            nkpts, nband = band.shape[-2:]
            band = band.reshape(-1, nkpts, nband) * const.AUTOEV
        
        # Determine efermi
        nkpts, nband = band.shape[-2:]
        nelectron = nkpts * self.nelec
        nocc = nelectron // 2       #TODO: uhf is considered 
        energy = band.flatten()
        idx = energy.argsort()
        efermi = energy[idx][:nocc].max()
        
        return efermi

    
    def get_band(self, band=None, kpts=None):
        '''Processing band structure info'''
        if band is None: 
            assert self.band is not None, "You need to provide bands calculated by PySCF kks.get_bands function"
            band = np.asarray(self.band[0])
            nkpts, nband = band.shape[-2:]
            band = band.reshape(-1, nkpts, nband) * const.AUTOEV
            
        if kpts is None: 
            assert self.kpts is not None, "You need to provide kpts used in the PySCF kks.get_bands function"
            kpath_frac = self.kpts     
        
        # Find absolute kpts and shift the band
        lattice = self.cell[0]
        recip_lattice = 2 * np.pi * np.linalg.inv(lattice).T # Get the reciprocal lattice in the row vector format
        abs_kpts = kpath_frac.dot(recip_lattice)                  # From fractional to absolute
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            
        # Get efermi
        efermi = self.get_efermi(band)
        
        return band, kpath_frac, proj_kpath, recip_lattice, efermi
        
    def get_pband(self, band=None, kpts=None):
        if band is None: 
            assert self.band is not None, "You need to provide bands calculated by PySCF kks.get_bands function"
            mo_coeff = np.asarray(self.band[1])
            nkpts, nband = band.shape[-2:]
            band = band.reshape(-1, nkpts, nband) * const.AUTOEV
            
        if kpts is None: 
            assert self.kpts is not None, "You need to provide kpts used in the PySCF kks.get_bands function"
            kpts = self.kpts     
            
        kpath_frac = kpts
        
        # Find absolute kpts and shift the band
        lattice = self.cell[0]
        recip_lattice = 2 * np.pi * np.linalg.inv(lattice).T # Get the reciprocal lattice in the row vector format
        abs_kpts = kpts.dot(recip_lattice)                  # From fractional to absolute
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            
        # Determine efermi
        nelectron = nkpts * self.nelec
        nocc = nelectron // 2       #TODO: uhf is considered 
        energy = band.flatten()
        idx = energy.argsort()
        efermi = energy[idx][:nocc].max()
        
        return band, kpath_frac, proj_kpath, recip_lattice, efermi
    
    def _generate_band(self, efermi=None, spin=0, klabel=None):
        '''Processing/collecting the band data before the plotting function
           klabel            : a list of labels and corresponding coordinates for high symmetry k-points
        '''          
        band, kpath_frac, proj_kpath, recip_lattice, efermi_ = self.get_band()
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
        # Find absolute coordinates for high symmetric kpoints  
        sym_kpoint_coor = None
        if klabel is not None:
            klabel, coor_kpts = str_format.format_klabel(klabel)  
            assert coor_kpts is not None, "You need to provide the coordinates for high symmetric k in the klabel"  
            abs_kpts = coor_kpts.dot(recip_lattice)   
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
                
        return band, proj_kpath, sym_kpoint_coor, klabel

    def get_bandgap(self):
        '''Get the bandgap'''       
        band, kpath_frac, proj_kpath, recip_lattice, efermi = self.get_band()
        nspin, nkpts, nbands = band.shape
        for spin in range(nspin):
            print('Spin:', spin)  
            CBM = None
            for bandth in range(nbands):
                shifted_band = band[spin,:,bandth] - efermi
                if (shifted_band > 0.0).all() == True:
                    CBM = band[spin,:, bandth]
                    VBM = band[spin,:, bandth -1]                
                    break
                elif ((shifted_band < 0.0).any() == True) and ((shifted_band > 0.0).any() == True):
                    print("This is a metal")
                    break
                    
            if CBM is not None:
                vbm_idx = np.argmax(VBM)
                cbm_idx = np.argmin(CBM)
                bandgap = CBM[cbm_idx] - VBM[vbm_idx]
                direct = False
                if vbm_idx == cbm_idx: direct = True

                print('  E(VBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (VBM[vbm_idx], 
                                                                kpath_frac[vbm_idx,0], kpath_frac[vbm_idx,1], kpath_frac[vbm_idx,2]))
                print('  E(CBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (CBM[cbm_idx], 
                                                                kpath_frac[cbm_idx,0], kpath_frac[cbm_idx,1], kpath_frac[cbm_idx,2]))
                if direct == True: 
                    print('  Direct bandgap   : %6.3f' % (bandgap))             
                else:  
                    print('  Indirect bandgap : %6.3f' % (bandgap))              
                    gap1 = CBM[vbm_idx] - VBM[vbm_idx]
                    gap2 = CBM[cbm_idx] - VBM[cbm_idx]
                    direct_gap = min(gap1, gap2)
                    print('  Direct bandgap   : %6.3f' % (direct_gap))
 
    def _generate_pband(self, band=None, spin=0, gradient=False, lm='spd'):
        '''Processing/collecting the projected band data before the plotting function
            
            Note: 
            proj_wf = [kpts, band, # of orbitals]
                  
            Examples for lm:
                lm = 'Ni:s ; p ; d'             :   three groups: (1) s of Ni ; (2) all p orbitals ; (3) all d orbitals
                lm = ['Ni:s,pd', 'O1:p;O2']     :   two groups: (1) s,p,d of Ni ; (2) all p orbitals of the 1st O  and all otbitals of O2 
                lm = ['Ni1;O', 'N']             :   two groups: (1) the 1st Ni and all the O atoms ; (2) All N atom
 
            if gradient == True: user has to provide a TWO groups of orbitals  
                for example, lm = 'Ni:s ; p' or ['Ni:s,pd', 'O1:p;O2']  
            
        '''      

        if band is None: 
            assert self.band is not None, "You need to provide bands calculated by PySCF kks.get_bands function"
            mo_coeff = np.asarray(self.band[1])         
            nkpts = mo_coeff.shape[-3]
            mo_coeff = mo_coeff.reshape(-1, nkpts, self.nao, self.nao) # dimension [spin, kpts, ao_th, band]
            mo_coeff = mo_coeff.transpose(0,1,3,2)  # dimension [spin, kpts, band, ao_th]
        
        proj_wf = np.absolute(mo_coeff[spin])
        species, lm_list =  pscf_utils.make_basis_dict(self.pyscf_cell)
        

        # Generate pband
        formatted_atom, formatted_lm = str_format.general_lm(lm)
        if gradient:        
            assert len(formatted_atom) == 2, "For the gradient plot, you only need to provide two groups of orbitals, for example, lm = 's,p'"

        # Calculate total band and pband
        total = proj_wf.sum(axis=2)     # sum over all orbitals, dimension now is [kpts, band]
        pband = []             
        for i, atoms in enumerate(formatted_atom):  
            proj_val = 0
            for j, atom in enumerate(atoms):
                # Locate the atom
                if atom is None: 
                    idx_atom = np.arange(len(species))
                else:
                    atom_, id = str_format.format_atom(atom)
                    assert atom_ in self.element, "This is wrong: " + atom + ". Check the lm string/list. Atom is must be in the element list: " + " ".join(self.element)
                    available_atom = [(n, atm) for n, atm in enumerate(self.atom) if atm == atom_]
                    natom = len(available_atom)
                    if id is not None:
                        assert id <= natom, "This is wrong: " + atom + ". Check the lm string/list. Atom id is must be <= " + str(natom) + " for: " + atom_
                        
                    idx_atom = []
                    nspecies = species.count(atom_)
                    nwfc = nspecies // natom
                    count = 0
                    for n, atm in enumerate(species):
                        if atm == atom_: 
                            if id is None:
                                idx_atom.append(n)
                            elif count // nwfc == id - 1: 
                                idx_atom.append(n)
                            count += 1 
                        
                # Locate the lm 
                idx_lm = []
                for idx in idx_atom:
                    for each_lm in formatted_lm[i][j]:
                        if each_lm is None:
                            idx_lm.append(idx)
                        elif lm_list[idx] == each_lm:
                            idx_lm.append(idx)              

                proj_val += (proj_wf[:,:,idx_lm]).sum(axis=2)
            pband.append(proj_val/total)
        pband = np.asarray(pband)
        
        if gradient:  
            pband = pband[0]/(pband.sum(axis=0))

        return pband   
