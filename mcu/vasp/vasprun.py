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
from ..cell import spg_wrapper, cell_io, cell
from ..cell import utils as cell_utils
from . import utils, vasp_io
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

        
class main(cell.main, plot.main):
    def __init__(self, prefix=None):
        '''
            prefix      : the working directory
        '''
        
        # Create vasprun object(s)
        if prefix is None: 
            self.prefix = os.getcwd()
        else:
            self.prefix = prefix
            
        if isinstance(self.prefix, str):                   # For one vasprun.xml file    
            self.vasprun = vasp_io.XML(self.prefix + '/vasprun.xml')
            self.get_info(self.vasprun)
            
        elif isinstance(self.prefix, list):                # For multiple vasprun.xml file
            self.vasprun = []
            for xml in self.prefix:
                xml_file = xml + '/vasprun.xml'
                assert check_exist(xml_file), 'Cannot find the vasprun.xml file. Check the prefix:' + xml_file
                self.vasprun.append(vasp_io.XML(xml_file))                
            self.get_info(self.vasprun[0])      # Only get info for the first vasprun.xml
        else:
            print('Provide a string or a list of names for *.xml file')
            
    def get_cell(self, vasprun):
        '''Get the cell info from vasprun and return the cell in spglib format'''
        self.cell_init  =  utils.cell_to_spgcell(vasprun.cell_init, self.atom)
        self.cell  = utils.cell_to_spgcell(vasprun.cell_final, self.atom)

    def get_info(self, vasprun):    
        '''Extract basis information from the vasprun.xml'''

        electronic = vasprun.parameters['electronic']
        self.nelec = electronic.general['NELECT']    
        self.nbands = electronic.general['NBANDS']
        self.soc = electronic.spin['LSORBIT']
        self.ispin = electronic.spin['ISPIN']    
        self.isym = vasprun.parameters['symmetry']['ISYM']
        self.isym_symprec = vasprun.parameters['symmetry']['SYMPREC']
        grids = vasprun.parameters['grids']
        self.ngrid = [grids['NGX'], grids['NGY'], grids['NGZ']]
        self.kmesh_type = vasprun.kpoints['type']
        if self.kmesh_type != 0 : self.kmesh = vasprun.kpoints['divisions']
        if self.kmesh_type == 2 :
            self.kmesh_shift = vasprun.kpoints['usershift']
        self.kpts = vasprun.kpoints['kpointlist']
        self.kpts_weight = vasprun.kpoints['weights']
        self.nkpts = self.kpts.shape[0] 
        self.natom  = vasprun.natom 
        self.atom  = vasprun.atom     
        self.element = [atom[1] for atom in vasprun.types]
        self.get_cell(vasprun)
        self.get_efermi()      
        
    def get_kpts_nosym(self, symprec=1e-5, no_spatial=False):    
        '''K-point list in VASP is an irreducible list
           This function maps the irreducible list to the full list using spglib
        '''
        
        assert self.kmesh_type == 'Monkhorst-Pack', \
                        "This function can only use for a wave function used a k-mesh scheme like Gamma center or Monkhorst-Pack"
        if self.isym == -1: 
            is_time_reversal = False
        else:
            is_time_reversal = True

        mapping_kpts, kpts_grid = self.get_mapping_kpts(mesh=self.kmesh, cell=self.cell, is_shift=self.kmesh_shift, symprec=symprec, 
                    is_time_reversal=is_time_reversal, no_spatial=no_spatial)

        return mapping_kpts, kpts_grid


############ Plotting #################
    def get_efermi(self):
        '''Extract E_fermi either from vasprun.xml or OUTCAR'''
        if isinstance(self.vasprun, vasp_io.XML):
            tdos, pdos, efermi = self.vasprun.get_dos()
            if efermi is not None:
                return efermi
            else:           # Get E_fermi from OUTCAR
                efermi = vasp_io.get_efermi_from_OUTCAR(self.prefix + "/OUTCAR")
                if efermi is not None:
                    return efermi
                else:
                    print("E_fermi cannot be read from vasprun.xml or OUTCAR, hence it is set to zero")
                    return 0
        elif isinstance(self.vasprun, list):        
            efermi_list = []
            for i in range(len(self.prefix)):
                tdos, pdos, efermi = self.vasprun[i].get_dos()
                if efermi is not None:
                    efermi_list.append(efermi)
                else:           # Get E_fermi from OUTCAR
                    efermi = vasp_io.get_efermi_from_OUTCAR(self.prefix + "/OUTCAR")
                    if efermi is None:
                        efermi = 0.0
                        print("E_fermi cannot be read from vasprun.xml or OUTCAR, hence it is set to zero")
                    efermi_list.append(efermi)
            return efermi_list
            
    def get_bandgap(self, efermi=None):
        '''Get the bandgap'''
        
        # Get the fermi level
        if efermi is None: efermi = self.get_efermi()
            
        if isinstance(self.vasprun, vasp_io.XML):              # For one vasprun.xml file
            band_occ = self.vasprun.get_band()
            band, co_occ = band_occ[:,:,:,0], band_occ[:,:,:,1]
            co_occ_ = co_occ > 0.5       
            electronic = self.vasprun.parameters['electronic']
        elif isinstance(self.vasprun,list):                             # For multiple vasprun.xml file
            assert isinstance(efermi,list), "There are more than one vasprun.xml, hence please provide a list of efermi"
            electronic = self.vasprun[0].parameters['electronic']
            nbands = electronic.general['NBANDS']
            
            band_spin = [] 
            co_occ_spin1 = []
            co_occ_spin2 = []            
            for spin in range(self.ispin):
                bands = np.zeros([1,nbands])
                co_occ1 = np.zeros([1,nbands])    
                co_occ2 = np.zeros([1,nbands], dtype=bool)                 
                kptss = np.zeros([1,3])
                for i, vasprun in enumerate(self.vasprun):          # loop over vasprun.xml
                    band_occ = vasprun.get_band()
                    band = band_occ[spin,:,:,0]
                    kpts = vasprun.kpoints['kpointlist']
                    weight = vasprun.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    kpts, band = kpts[nonzero:], band[nonzero:]
                    co_occ = band_occ[spin,nonzero:,:,1]
                    co_occ_ = band < efermi[i] 
                    bands = np.vstack([bands,band])
                    kptss = np.vstack([kptss,kpts])
                    co_occ1 = np.vstack([co_occ1,co_occ])
                    co_occ2 = np.vstack([co_occ2,co_occ_])
                band_spin.append(bands[1:])  
                co_occ_spin1.append(co_occ1[1:])   
                co_occ_spin2.append(co_occ2[1:])  
            self.kpts, band = np.asarray(kptss[1:]), np.asarray(band_spin)
            self.nkpts = self.kpts.shape[0]
            co_occ, co_occ_ = np.asarray(co_occ_spin1), np.asarray(co_occ_spin2)
            
        bandedge = np.zeros([self.ispin,self.nkpts,2,2])
        self.bandgap = []
        for spin in range(self.ispin):
            print('Spin:', spin)        
            for kpt in range(self.nkpts):
                band_kpt = band[spin,kpt]
                occ = co_occ_[spin,kpt]               
                homo_idx = np.count_nonzero(occ) - 1
                lumo_idx = homo_idx + 1               
                bandedge[spin,kpt,0,0] = band_kpt[homo_idx]
                bandedge[spin,kpt,0,1] = co_occ[spin,kpt,homo_idx]
                bandedge[spin,kpt,1,0] = band_kpt[lumo_idx]
                bandedge[spin,kpt,1,1] = co_occ[spin,kpt,lumo_idx]
                
            vbm_idx = np.argmax(bandedge[spin,:,0,0])
            cbm_idx = np.argmin(bandedge[spin,:,1,0])
            direct = False
            if vbm_idx == cbm_idx: direct = True
            print('  E(VBM) = %7.4f with occ = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (bandedge[spin,vbm_idx,0,0], bandedge[spin,vbm_idx,0,1], 
                                                            self.kpts[vbm_idx,0], self.kpts[vbm_idx,1], self.kpts[vbm_idx,2]))
            print('  E(CBM) = %7.4f with occ = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (bandedge[spin,cbm_idx,1,0], bandedge[spin,cbm_idx,1,1], 
                                                            self.kpts[cbm_idx,0], self.kpts[cbm_idx,1], self.kpts[cbm_idx,2]))
            bandgap = bandedge[spin,cbm_idx,1,0] - bandedge[spin,vbm_idx,0,0] 
            self.bandgap.append(bandgap)  
            if direct == True: 
                print('  Direct bandgap   : %6.3f' % (bandgap))             
            else:  
                print('  Indirect bandgap : %6.3f' % (bandgap))              
                gap1 = bandedge[spin,cbm_idx,1,0] - bandedge[spin,cbm_idx,0,0]
                gap2 = bandedge[spin,vbm_idx,1,0] - bandedge[spin,vbm_idx,0,0]  
                direct_gap = gap1
                direct_k = self.kpts[cbm_idx]
                if gap1 > gap2: 
                    direct_gap = gap2
                    direct_k = self.kpts[vbm_idx]
                    
                print('  Direct bandgap   : %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (direct_gap, direct_k[0], direct_k[1], direct_k[2]))                   
                
                
    def _generate_band(self, vasprun=None, efermi=None, spin=0, klabel=None):
        '''Processing/collecting the band data before the plotting function
        '''
        if vasprun is None: vasprun = self.vasprun
        
        # Get the fermi level
        if efermi is None: efermi = self.get_efermi()   
        
        sym_kpoint_coor = None
        band = None
        proj_kpath = None            #projected kpath: kpath coordinated projected to 1D
            
        if isinstance(vasprun, vasp_io.XML) and vasprun.kpoints['type'] == 1: # For conventional band structure calculation 
            if klabel is not None:
                klabel, coor_kpts = str_format.format_klabel(klabel)
            
            band = vasprun.get_band()[spin,:,:,0]
            kpts = vasprun.kpoints['kpointlist']
            kpts, band = utils.rm_redundant_band(kpts, band) 
            
            # Find absolute kpts and shift the band
            b = vasprun.cell_final[1]               # Get the reciprocal lattice in the row vector format
            abs_kpts = kpts.dot(b)                  # From fractional to absolute
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            band = band - efermi               # Efermi is set at 0 eV
            
            highsym_kpt = vasprun.kpoints['points']
            nkpts = highsym_kpt.shape[0]
            sym_kpoint_coor = [0.0]
            for kpt in range(nkpts-2):
                idx = ((proj_kpath.shape[1] + nkpts - 2)//(nkpts-1) - 1) * (kpt+1)
                coor = proj_kpath[0,idx]         
                sym_kpoint_coor.append(coor)
            sym_kpoint_coor.append(1.0*proj_kpath.max())   
            sym_kpoint_coor = np.asarray(sym_kpoint_coor)
            if klabel is not None:
                assert len(klabel) == len(sym_kpoint_coor), "The number of k label must be " + str(len(sym_kpoint_coor))
        else:
            if isinstance(vasprun, vasp_io.XML):                       # For one vasprun.xml file
                band = vasprun.get_band()[spin,:,:,0]
                kpts = vasprun.kpoints['kpointlist']
                if vasprun.kpoints['type'] == 0:
                    weight = vasprun.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    kpts, band = kpts[nonzero:], band[nonzero:]
                band = band - efermi
            elif isinstance(vasprun,list):                                      # For multiple vasprun.xml file
                assert isinstance(efermi,list)
                electronic = vasprun[0].parameters['electronic']
                nbands = electronic.general['NBANDS']
                bands = np.zeros([1,nbands])
                kptss = np.zeros([1,3])
                for i, run in enumerate(vasprun):
                    band = run.get_band()[spin,:,:,0]
                    kpts = run.kpoints['kpointlist']
                    weight = run.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    kpts, band = kpts[nonzero:], band[nonzero:]
                    band = band - efermi[i]
                    bands = np.vstack([bands,band])
                    kptss = np.vstack([kptss,kpts])
                    
                kpts, band = kptss[1:], bands[1:]
                vasprun = vasprun[0]
                
            # Find absolute kpts
            b = vasprun.cell_final[1]               # Get the reciprocal lattice
            abs_kpts = kpts.dot(b)                  # From fractional to absolute in A^-1 unit
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())

            # Find absolute coordinates for high symmetric kpoints  
            if klabel is not None:
                klabel, coor_kpts = str_format.format_klabel(klabel)  
                assert coor_kpts is not None, "You need to provide the coordinates for high symmetric k in the klabel"    
                abs_kpts = coor_kpts.dot(b)   
                temp_kpts = np.empty_like(abs_kpts)
                temp_kpts[0] = abs_kpts[0]
                temp_kpts[1:] = abs_kpts[:-1] 
                sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
        
        return band, proj_kpath, sym_kpoint_coor, klabel
                
    def _generate_pband(self, vasprun=None, spin=0, gradient=False, lm='spd'):
        '''Processing/collecting the projected band data before the plotting function
            proj_wf = [kpt,band,atom,lm] , read vasp_io.XML.get_projected for more details info
            
            Examples for lm:
                lm = 'Ni:s ; p ; d'             :   three groups: (1) s of Ni ; (2) all p orbitals ; (3) all d orbitals
                lm = ['Ni:s,pd', 'O1:p;O2']     :   two groups: (1) s,p,d of Ni ; (2) all p orbitals of the 1st O  and all otbitals of O2 
                lm = ['Ni1;O', 'N']             :   two groups: (1) the 1st Ni and all the O atoms ; (2) All N atom
 
            if gradient == True: user has to provide a TWO groups of orbitals  
                for example, lm = 'Ni:s ; p' or ['Ni:s,pd', 'O1:p;O2']               
        '''      
        if vasprun is None: vasprun = self.vasprun
        
        # Collecting/combining the projected wfn from vasprun.xml
        if isinstance(vasprun, vasp_io.XML):                       # For one vasprun.xml file
            vasprun.get_projected()
            proj_wf = vasprun.proj_wf[spin] 
            lm_list = vasprun.lm           
            proj_wf = utils.rm_redundant_band(self.kpts, proj_wf)[1]          # remove redundant 
        
        elif isinstance(vasprun,list):                                      # For multiple vasprun.xml file
            nlm    = len(vasprun[0].get_lm())
            lm_list = vasprun[0].lm
            proj_wfs = np.zeros([1,self.nbands,self.natom,nlm])
            
            for i, run in enumerate(vasprun):
                run.get_projected()
                proj_wf = run.proj_wf[spin]
                kpts = run.kpoints['kpointlist']
                weight = run.kpoints['weights']
                nonzero = np.count_nonzero(weight)
                proj_wf = proj_wf[nonzero:]
                proj_wfs = np.concatenate((proj_wfs,proj_wf))   
            proj_wf = proj_wfs[1:]
         
        # Generate pband
        formatted_atom, formatted_lm = str_format.general_lm(lm)
        
        if gradient:        
            assert len(formatted_atom) == 2, "For the gradient plot, you only need to provide two groups of orbitals, for example, lm = 's,p'"

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
                
    def _generate_dos(self, vasprun=None, efermi=None, spin=0, lm=None):
        '''Processing/collecting the DOS data before the plotting function
            Note: unlike plot_band function, only one vasprun.xml is used. Combining vasprun.xml for DOS sounds a bit weird
            and unececessary
            
            pdos_data = [E(eV), atom-th, lm-th]
            spin            : spin of DOS.
                              For LSORBIT == True: spin = 0,1,2,3
                              For ISPIN = 2      : spin = 0,1
                              
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        if lm is None: 
            lm = [atom+':s,p,d' for atom in self.element] 
            
        if vasprun is None: 
            if isinstance(self.vasprun, vasp_io.XML): 
                vasprun = self.vasprun
                if efermi is None: efermi = self.get_efermi() 
            if isinstance(self.vasprun,list): 
                vasprun = self.vasprun[0]  
                if efermi is None: efermi = self.get_efermi()[0]
        else:
            assert isinstance(vasprun, vasp_io.XML)
            
        tdos, pdos, e_fermi = vasprun.get_dos()
        assert tdos is not None, "Cannot find DOS data in vasprun.xml" 
        tdos = tdos[spin,:,:2]
        lm_list = vasprun.lm
        
        # Compute pDOS
        if pdos is not None:
            pdos_spin = pdos[spin,:,:,1:]
            formatted_atom, formatted_lm = str_format.general_lm(lm)

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
                    
                    proj_val += (pdos_spin[:,idx_atom, :][:,:,idx_lm]).sum(axis=(1,2))
                pdos.append(proj_val)
            pdos = np.asarray(pdos).T 
  
        # Shift the energy 
        tdos[:,0] = tdos[:,0] - efermi

        return tdos, pdos
        
    def _generate_spin(self, vasprun, lm=None):
        '''Processing/collecting the spin texture data before the plotting function
        
            proj_wf = [spin, kpt, band, atom, lm] 
            read vasp_io.XML.get_projected for more details info
            
            lm contains a group of orbitals
            Examples for lm:
                lm = 'spd'
                lm = [' Ni:s,pd ; O1:p ; O2 ']
        '''      
        
        formatted_atom, formatted_lm = str_format.general_lm(lm)
        assert len(formatted_atom) == 1, "For spin texture plot, you only need to provide one group of orbitals, for example, lm = 's,p'"

        # Collecting/combining the projected wfn from vasprun.xml
        if isinstance(vasprun, vasp_io.XML):                       # For one vasprun.xml file
            vasprun.get_projected()
            lm_list = vasprun.lm 
            proj_wfs = []
            for spin in range(1,4):
                proj_wf = vasprun.proj_wf[spin] 
                nonzero = np.count_nonzero(self.kpts_weight)
                proj_wf = proj_wf[nonzero:]
                proj_wfs.append(proj_wf)
                
        elif isinstance(vasprun,list):                                      # For multiple vasprun.xml file
            nlm    = len(vasprun[0].get_lm())
            lm_list = vasprun[0].lm
            proj_wfs = []
            for spin in range(1,4):
                proj_wf = np.zeros([1,self.nbands,self.natom,nlm])         
                for i, run in enumerate(vasprun):
                    run.get_projected()
                    proj = run.proj_wf[spin]
                    weight = run.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    proj = proj[nonzero:]
                    proj_wf = np.concatenate((proj_wf,proj))  
                proj_wf = proj_wf[1:]
                proj_wfs.append(proj_wf)
            
        spin_text = []  
        for spin in range(3):        
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
                    
                    proj_val += (proj_wfs[spin][:,:,idx_atom, :][:,:,:,idx_lm]).sum(axis=(2,3))
            spin_text.append(proj_val)
        spin_text = np.asarray(spin_text)
        
        return spin_text

    def plot_spin(self, style=1, lm=None, band=None, cmap='bwr', color='k', scale=15, scale_units=None,
                    save=False, figname='spin_texture', figsize=(7,6), xlim=None, ylim=None, fontsize=18, dpi=600, format='png'):
        '''Plot spin texture
           
            Attribute:
                efermi      : a Fermi level or a list of Fermi levels, it is automatically extracted frim vasprun.xml or OUTCAR
                lm          : 'Ni:s' or ['Ni:s','C:pz']
                
                band        : index of the band
                color       : color of the arrow for style 1
                scale       : the size of the marker
                alpha       : the transparency level of curves
                cmap        : color map in the style 3, following the matplotlib
                edgecolor   : the marker's border color in the style 3, default: 'none', any color code should work
                facecolor   : the filling color of style 1 and 2
                              None              : taking from the color arg
                              'none'            : unfilling circle
                              [False,True,...]  : True for filling markers and False for empty ones
                marker      : a list of marker shape, default is: 'o'
                legend      : a list of labels for different group of orbitals (same color) for the style 1 and 2            
        '''

        if lm == None : lm = [atom+':s,p,d' for atom in self.element]         
        
        # Check if the band values are reasonable otherwise generate it
        if band == None:
            idx_vbm = int(self.nelec)
            if self.soc == False: idx_vbm = idx_vbm//2
            band = idx_vbm - 1
        else:
            assert (band <= self.nbands) and (band >= 1)  
            band = band - 1

        spin_text = self._generate_spin(self.vasprun, lm=lm)[:,:,band]

        # Get X, Y
        kpoint = vasp_io.KPOINTS()
        plane, krange, npoint = kpoint.get_spin_kmesh()
        deltax = 2*krange[0]/(npoint[0]-1)
        deltay = 2*krange[1]/(npoint[1]-1)
        X, Y = np.mgrid[-krange[0]:(krange[0]+deltax*0.5):deltax,-krange[1]:(krange[1]+deltay*0.5):deltay]
        spin_text = spin_text.reshape(3,npoint[0],npoint[1])
        self.spin_text = spin_text
        if plane == 'xy':
            U = spin_text[0]
            V = spin_text[1]
            C = spin_text[2]
        elif plane == 'xz':
            U = spin_text[0]
            V = spin_text[2]
            C = spin_text[1]
        elif plane == 'yz':        
            U = spin_text[1]
            V = spin_text[2]
            C = spin_text[0]  
        
        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        if xlim == None: xlim = [X.min(),X.max()]
        if ylim == None: ylim = [Y.min(),Y.max()]
        yx_ratio = (ylim[1]-ylim[0]) /(xlim[1]-xlim[0])         
            
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        if style == 1:
            plt.contour(X, Y, C, linewidths=0.01)
            plt.contourf(X, Y, C, vmin=-1.0, vmax=1.0, cmap=cmap, levels=np.linspace(-1.0,1.0,1000))
            cbar = plt.colorbar(ticks=[-1.0,0.0,1.0], fraction=0.047*yx_ratio, pad=0.04*yx_ratio)  # Magic number that adjust the colorbar size        
            plt.quiver(X, Y, U, V, color=color, scale=scale, scale_units=scale_units)
        elif style == 2:
            plt.quiver(X, Y, U, V, C, cmap=cmap, scale=scale, scale_units=scale_units)
            cbar = plt.colorbar(ticks=[-1.0,0.0,1.0], fraction=0.047*yx_ratio, pad=0.04*yx_ratio)   # Magic number that adjust the colorbar size
           
        # Graph adjustments  
        border = 1.08
        ax.tick_params(labelsize=fontsize, width=border)
        ax.spines['top'].set_linewidth(border)
        ax.spines['right'].set_linewidth(border)
        ax.spines['bottom'].set_linewidth(border)
        ax.spines['left'].set_linewidth(border)
        ax.tick_params(labelsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.clim(-1,1)        
        if plane == 'xy': 
            x_label = r'$k_x$'
            y_label = r'$k_y$'   
        elif plane == 'xz':
            x_label = r'$k_x$'
            y_label = r'$k_z$'              
        elif plane == 'yz':
            x_label = r'$k_y$'
            y_label = r'$k_z$'        
        ax.set_xlabel(x_label + r' ($\AA^{-1}$)', size=fontsize+4)
        ax.set_ylabel(y_label + r' ($\AA^{-1}$)', size=fontsize+4)      
        ax.set_xticks([xlim[0],0,xlim[1]])
        ax.set_yticks([ylim[0],0,ylim[1]])
        plt.gca().set_aspect('equal', adjustable='box')
        plt.tight_layout()
        if save == True: 
            fig.savefig(figname+'.'+format,dpi=dpi,format=format)      
        else:
            plt.show()   
                
    def plot_band2D(self, efermi=None, spin=0, band=None, cmap='jet', save=False,
                    figsize=(8,6), figname='BAND2D', xlim=None, ylim=None, zlim=None, fontsize=16, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                                  0 or 1 for spin polarized
                label           : label for high symmetric points, e.g. 'G-X-L-W-G'
                                  if hybridXC=True, the lavel should be a list of labels plus their coordinates
                color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                                  
                                  
        '''
        
        band_idx = band
        if band_idx == None:
            idx_vbm = int(self.nelec)
            if self.soc == False: idx_vbm = idx_vbm//2               # Estimation for insulator, generally wrong for metal
            band_idx = [idx_vbm-1, idx_vbm]
        else:
            assert band_idx[0] <= band_idx[1]                    # from band[0] to band[1]
            if band_idx[0] < 1: band_idx[0] = 1     # index used in OUTCAR, will be shifted to start at zero
            if band_idx[1] > self.nbands: band_idx[1] = self.nbands              # Cannot larger than the number of bands
            band_idx[0] = band_idx[0] -1
            band_idx[1] = band_idx[1] -1            
            
        band, proj_kpath, sym_kpoint_coor, klabel = self._generate_band(self.vasprun, efermi, spin, klabel=None)  
        
        # Get X, Y
        kpoint = vasp_io.KPOINTS()
        plane, krange, npoint = kpoint.get_spin_kmesh()
        deltax = 2*krange[0]/(npoint[0]-1)
        deltay = 2*krange[1]/(npoint[1]-1)
        X, Y = np.mgrid[-krange[0]:(krange[0]+deltax*0.5):deltax,-krange[1]:(krange[1]+deltay*0.5):deltay]
        Z = band.reshape(npoint[0],npoint[1],self.nbands)
            
        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')

        # Plot bands   
        vmin = Z[:,:,band_idx[0]:band_idx[1]+1].min()
        vmax = Z[:,:,band_idx[0]:band_idx[1]+1].max()
        for ith in range(band_idx[0],band_idx[1]+1):
            surf = ax.plot_surface(X, Y, Z[:,:,ith], cmap=cmap, vmin=vmin, vmax=vmax, linewidth=0, antialiased=False) 
        
        # Graph adjustments             
        ax.tick_params(labelsize=fontsize)
        if xlim == None: xlim = [X.min(),X.max()]
        if ylim == None: ylim = [Y.min(),Y.max()]
        if zlim != None: plt.zlim(zlim)
        plt.xlim(xlim)
        plt.ylim(ylim)
        if plane == 'xy': 
            x_label = r'$k_x$'
            y_label = r'$k_y$'   
        elif plane == 'xz':
            x_label = r'$k_x$'
            y_label = r'$k_z$'              
        elif plane == 'yz':
            x_label = r'$k_y$'
            y_label = r'$k_z$'        
        ax.set_xlabel(x_label + r' ($\AA^{-1}$)', size=fontsize+4)
        ax.set_ylabel(y_label + r' ($\AA^{-1}$)', size=fontsize+4)
        ax.set_zlabel('Energy (eV)', size=fontsize+4)        
        ax.set_xticks([xlim[0],0,xlim[1]])
        ax.set_yticks([ylim[0],0,ylim[1]])
        ax.xaxis.labelpad = 15
        ax.yaxis.labelpad = 15
        ax.zaxis.labelpad = 15
        cbar = fig.colorbar(surf, shrink=0.5, aspect=20)      
        cbar.ax.tick_params(labelsize=fontsize+4)
        plt.tight_layout()
        if save == True: 
            fig.savefig('Band.'+format,dpi=dpi,format=format)      
        else:
            plt.show() 
            