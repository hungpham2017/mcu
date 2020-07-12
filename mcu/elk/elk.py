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
from ..utils import plot, str_format
from ..utils.misc import check_exist
from ..vasp import const
from ..cell import spg_wrapper, cell_io, cell
from ..cell import utils as cell_utils
from . import elk_io
        
class main(cell.main, plot.main):
    def __init__(self,  prefix="./"):
        '''
           Args:
            prefix          : the working directory path
        '''
        assert prefix is not None, "Provide a prefix name for your project, for example, prefix.out, prefix.bs, etc."
        self.prefix = prefix
        self.get_info() 

############ General #################
        
    def get_info(self, prefix=None):    
        '''Extract basis information from the vasprun.xml'''
        if prefix is None: prefix = self.prefix
        
        assert check_exist(self.prefix + "/INFO.OUT"), "Cannot find any INFO.OUT file"
        assert check_exist(self.prefix + "/KPOINTS.OUT"), "Cannot find any KPOINTS.OUT file"
        assert check_exist(self.prefix + "/EIGVAL.OUT"), "Cannot find any KPOINTS.OUT file"
        
        data = elk_io.read_info(self.prefix + "/INFO.OUT")
        self.cell = data['cell']
        self.atom = data['atom']
        self.natom  = len(self.atom)
        self.element = [atm['atom type'] for atm in data['Species']]
        self.nelec = int(data['charge']['valence'])
        self.nbands = int(data['states']['valence'])
        self.ispin = data['spin'] 
        self.soc = data['soc']
        self.efermi = data['Energies'][-1]['Fermi']
        self.kpts, self.kpts_weight = elk_io.read_kpoints(self.prefix + "/KPOINTS.OUT")
        self.eigvals, self.occ = elk_io.read_eigval(self.prefix + "/EIGVAL.OUT")
        
############ Plotting ################# 
    def get_efermi(self):
        '''Get Fermi energy'''
        return self.efermi
        
    def get_band(self, filename=None):
        ''' Get band from the BAND.OUT'''
        
        found = check_exist(self.prefix + "/BAND.OUT") or check_exist(self.prefix + "/BAND_S01_A0001.OUT")
        assert found, "Cannot find any BAND.OUT file"
        
        if check_exist(self.prefix + "/BAND.OUT"):
            filename = self.prefix + "/BAND.OUT"
            proj_kpath, band = elk_io.read_band(filename, spin_polarized=self.ispin)
        elif check_exist(self.prefix + "/BAND_S01_A0001.OUT"):
            # this for task 21, 22, 23
            filename = self.prefix + "/BAND_S01_A0001.OUT"
            proj_kpath, band = elk_io.read_band(filename, spin_polarized=self.ispin)
        
        band = band[:,:,:,0] + self.efermi           
        lattice = self.cell[0]
        recip_lattice = 2 * np.pi * np.linalg.inv(lattice).T
        efermi = self.get_efermi()

        # Get the fractional high symmetric kpts from input file
        assert check_exist(self.prefix + "/elk.in"), "Cannot find any elk.in file"
        sym_kvecs, npts = elk_io.read_plot1d(self.prefix + "/elk.in")
        nkvecs = sym_kvecs.shape[0]

        # Compute the fractional kpts
        abs_kpoint = sym_kvecs.dot(recip_lattice)
        temp_kpts = np.empty_like(abs_kpoint)
        temp_kpts[0] = abs_kpoint[0]
        temp_kpts[1:] = abs_kpoint[:-1] 
        ksegment = np.sqrt(((temp_kpts - abs_kpoint)**2).sum(axis=1))
        ksegment_cumsum = ksegment.cumsum()
        self.sym_kpoint_coor = ksegment_cumsum
        
        total_pts = npts - sym_kvecs.shape[0]
        dt = ksegment_cumsum[-1]
        kpts = [sym_kvecs[0]]
        for i in range(nkvecs - 1):
            segment_pts = int(round(total_pts * ksegment[i + 1] /dt))
            if (segment_pts >= total_pts) or (i == nkvecs - 2): segment_pts = total_pts    
            kvecs_diff = [sym_kvecs[i + 1] - sym_kvecs[i]]
            kvecs = np.asarray(kvecs_diff * (segment_pts + 1)) 
            kvecs *= np.asarray([(np.arange(segment_pts + 1) + 1 )/ (segment_pts + 1)] * 3).T
            kvecs +=sym_kvecs[i]
            kpts.append(kvecs)
            dt = dt - ksegment[i + 1]
            total_pts = total_pts - segment_pts
        
        # sanity check to make sure the calculated kpts is what ELK is used in the band calculations
        kpts = np.vstack(kpts)
        abs_kpts = kpts.dot(recip_lattice)
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        kpath = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum()
        
        assert abs(kpath - proj_kpath).max() < 1e-8, "The calculated kpts does not mat the one used by ELK"

        return band, kpts, proj_kpath, recip_lattice, efermi

    def get_bandgap(self, filename=None):
        '''Get the bandgap'''
        
        band, kpath_frac, proj_kpath, recip_lattice, efermi = self.get_band(filename) 
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
                    
    def _generate_band(self, filename=None, efermi=None, spin=0, klabel=None):
        '''Processing/collecting the band data before the plotting function
           klabel            : a list of labels and corresponding coordinates for high symmetry k-points
        '''
        band, kpath_frac, proj_kpath, recip_lattice, efermi_ = self.get_band(filename)
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
        # Find absolute coordinates for high symmetric kpoints  
        sym_kpoint_coor = self.sym_kpoint_coor
        if klabel is not None:
            klabel, coor_kpts = str_format.format_klabel(klabel)  
            
        return band, proj_kpath, sym_kpoint_coor, klabel
        
    def _generate_pband(self, prefix=None, spin=0, gradient=False, lm='spd'):
        ''' Processing/collecting the projected band data before the plotting function
            In ELK, there are type options to project the band structure
                  
                task 21   : 
                    the contribution from s, p, d, f  (sum over m)  
                    
                task 22   : 
                    the contribution from s, p(-1,0,1), d(-2,-1,0,+1,+2), f(-3,-2,-1,0,+1,+2,+3)  
                    There is no the x,y,z-like orbitals (e.g., p, py, pz) ... yet.
                    TODO: either implement it in ELK (easier option) 
                          or find a work-around solution in MCU
                    IN my current work-around solution, let's just assume:
                        For p: 
                            m = -1 , 0 , +1   
                            equivalent to: px, pz, py
                        For d: 
                            m = -2, -1 , 0 , +1, +2   
                            equivalent to: dx2-y2, dxz, dz2, dyz, dxy
                        For f: 
                            m = -3, -2, -1 , 0 , +1, +2, +3   
                            equivalent to: fx(x2−3y2), fxyz, fxz2, dz3, fyz2, fz(x2−y2), fy(3x2−y2)                         

                task 23   : 
                    the contribution from the spin up and down in the spinor wave function.
                    This is useful for the case of SOC calculation
                    In this case, the lm argument will be diregarded, 
                    
                  
            Note:
                proj_wf = [kpt,band,atom,lm]
                
            Examples for lm:
                lm = 'Ni:s ; p ; d'             :   three groups: (1) s of Ni ; (2) all p orbitals ; (3) all d orbitals
                lm = ['Ni:s,pd', 'O1:p;O2']     :   two groups: (1) s,p,d of Ni ; (2) all p orbitals of the 1st O  and all otbitals of O2 
                lm = ['Ni1;O', 'N']             :   two groups: (1) the 1st Ni and all the O atoms ; (2) All N atom
 
            if gradient == True: user has to provide a TWO groups of orbitals  
                for example, lm = 'Ni:s ; p' or ['Ni:s,pd', 'O1:p;O2']  
                
            Return:
                pband
            
        '''      
        if prefix is None: prefix = self.prefix

        # Collect the BAND_Sss_Aaaaa files
        pband_data = []
        for i, element in enumerate(self.element):
            for atm in range(self.natom):
                filename = prefix + "/BAND_" + "S%02d_A%04d.OUT" % (i + 1, atm + 1)
                found = check_exist(filename)
                if found:
                    proj_kpath, atm_pband_data = elk_io.read_band(filename, spin_polarized=self.ispin)
                    pband_data.append(atm_pband_data[spin,:,:,1:])
                else:
                    break
            else:
                continue
        
        proj_wf = np.asarray(pband_data).transpose(1, 2, 0, 3) # [atom, kpt,band,lm] --> [kpt,band,atom,lm]
        # Generate pband
        if proj_wf.shape[-1] == 16:         # tast 22
            lm_list = ['s', 'px', 'pz', 'py', 'dx2-y2', 'dxz', 'dz2', 'dyz', 'dxy', 'fx(x2−3y2)', 'fxyz', 'fxz2', 'dz3', 'fyz2', 'fz(x2−y2)', 'fy(3x2−y2)']
            formatted_atom, formatted_lm = str_format.general_lm(lm)
            if gradient:        
                assert len(formatted_atom) == 2, "For the gradient plot, you only need to provide two groups of orbitals, for example, lm = 's,p'"
        elif proj_wf.shape[-1] == 5:         # tast 21
            proj_wf = proj_wf[:,:,:,1:]          # remove the total l column
            lm_list = ['s', 'p', 'd', 'f']
            formatted_atom, formatted_lm = str_format.general_lm(lm, type='l')
            if gradient:        
                assert len(formatted_atom) == 2, "For the gradient plot, you only need to provide two groups of orbitals, for example, lm = 's,p'"
        elif proj_wf.shape[-1] == 2:         # tast 23  for spinnor wave function
            lm_list = ['u', 'd']
            formatted_atom, formatted_lm = str_format.general_lm(lm, type='ud')
        elif proj_wf.shape[-1] == 1:         # tast 23 for spin-restricted wave function
            lm_list = ['s']
            formatted_atom, formatted_lm = str_format.general_lm(lm, type='spin')

        # Calculate total band and remove zero values
        total = proj_wf.sum(axis=(2,3))     # sum over all atoms and orbitals
        shape = total.shape
        idx_zeros = total.flatten() < 0.0001
        total = total.flatten()
        total[idx_zeros] = 1.0
        total = total.reshape(shape)

        pband = []             
        for i, atoms in enumerate(formatted_atom):  
            proj_val = 0
            for j, atom in enumerate(atoms):
                # Locate the atom
                if atom is None: 
                    idx_atom = np.arange(len(self.atom))
                else:
                    atom_, id = str_format.format_atom(atom)
                    assert atom_ in self.element, "This is wrong: " + atom + ". Check the lm string/list. Atom is must be in the element list: " + " ".join(self.element)
                    available_atom = [(n, atm) for n, atm in enumerate(self.atom) if atm == atom_]
                    natom = len(available_atom)
                    if id is not None:
                        assert id <= natom, "This is wrong: " + atom + ". Check the lm string/list. Atom id is must be <= " + str(natom) + " for: " + atom_
                        
                    idx_atom = []
                    nspecies = self.atom.count(atom_)
                    nwfc = nspecies // natom
                    count = 0
                    for n, atm in enumerate(self.atom):
                        if atm == atom_: 
                            if id is None:
                                idx_atom.append(n)
                            elif count  == id - 1: 
                                idx_atom.append(n)
                            count += 1 
                        
                # Locate the lm 
                idx_lm = []
                for idx, lm in enumerate(lm_list):
                    for each_lm in formatted_lm[i][j]:
                        if each_lm is None:
                            idx_lm.append(idx)
                        elif each_lm == lm:
                            idx_lm.append(idx)              
                
                proj_val += (proj_wf[:,:,idx_atom, :][:,:,:,idx_lm]).sum(axis=(2,3))
            pband.append(proj_val/total)
        pband = np.asarray(pband)
        
        if gradient:  
            pband = pband[0]/(pband.sum(axis=0))

        return pband 
        
    def _generate_dos(self, prefix=None, efermi=None, spin=0, lm=None):
        '''Processing/collecting the DOS data before the plotting function
            
            TDOS dimensions: [spin , [E(eV), tdos(E)]]
            
            spin            : spin of DOS.
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        
        if prefix is None: prefix = self.prefix    
        if self.soc: 
            spin_polarized = True
        else:
            spin_polarized = self.ispin
        if efermi is None: 
            efermi = self.get_efermi()
        if lm is None: 
            lm = [atom+':s,p,d,f' for atom in self.element] 
            
        # Get total DOS
        tdos_file = prefix + "/TDOS.OUT" 
        assert check_exist(tdos_file), "Cannot find " + tdos_file
        epsilon, tdos_data = elk_io.read_dos(tdos_file, spin_polarized=spin_polarized)
        tdos = np.asarray([epsilon, tdos_data[spin,:,0]]).T
        
        # Collect the pdos files
        pdos_data = []
        for i, element in enumerate(self.element):
            for atm in range(self.natom):
                filename = prefix + "/PDOS_" + "S%02d_A%04d.OUT" % (i + 1, atm + 1)
                found = check_exist(filename)
                if found:
                    lm_pdos_data = elk_io.read_dos(filename, spin_polarized=spin_polarized)[1]
                    pdos_data.append(lm_pdos_data[spin])
                else:
                    break
            else:
                continue
        
        pdos_data = np.asarray(pdos_data).transpose(1,0,2) # [atom, e, dos] --> e, atom, dos]
        
        # Create the possible lm list
        lm_list = ['s', 'px', 'pz', 'py', 'dx2-y2', 'dxz', 'dz2', 'dyz', 'dxy', 'fx(x2−3y2)', 'fxyz', 'fxz2', 'dz3', 'fyz2', 'fz(x2−y2)', 'fy(3x2−y2)']
        formatted_atom, formatted_lm = str_format.general_lm(lm)

        # Compute pDOS
        pdos = []             
        for i, atoms in enumerate(formatted_atom):  
            proj_val = 0
            for j, atom in enumerate(atoms):
                # Locate the atom
                if atom is None: 
                    idx_atom = np.arange(len(self.atom))
                else:
                    atom_, id = str_format.format_atom(atom)
                    assert atom_ in self.element, "This is wrong: " + atom + ". Check the lm string/list. Atom is must be in the element list: " + " ".join(self.element)
                    available_atom = [(n, atm) for n, atm in enumerate(self.atom) if atm == atom_]
                    natom = len(available_atom)
                    if id is not None:
                        assert id <= natom, "This is wrong: " + atom + ". Check the lm string/list. Atom id is must be <= " + str(natom) + " for: " + atom_
                        
                    idx_atom = []
                    count = 0
                    for n, atm in enumerate(self.atom):
                        if atm == atom_: 
                            if id is None:
                                idx_atom.append(n)
                            elif count  == id - 1: 
                                idx_atom.append(n)
                            count += 1 
                        
                # Locate the lm 
                idx_lm = []
                for idx, lm in enumerate(lm_list):
                    for each_lm in formatted_lm[i][j]:
                        if each_lm is None:
                            idx_lm.append(idx)
                        elif each_lm == lm:
                            idx_lm.append(idx)              
                
                proj_val += (pdos_data[:,idx_atom, :][:,:,idx_lm]).sum(axis=(1,2))
            pdos.append(proj_val)
        pdos = np.asarray(pdos).T 

        # Shift the energy 
        tdos[:,0] = tdos[:,0] - efermi 

        return tdos, pdos
        