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
from ..vasp import const
from ..cell import utils as cell_utils
from ..cell import cell
from . import qe_io

        
class main(cell.main, plot.main):
    def __init__(self,  prefix=None):
        assert prefix is not None, "Provide a prefix name for your project, for example, prefix.scf.out, prefix.band.out, etc."
        self.prefix = prefix
        self.get_info() 
        
############ General #################
        
    def get_info(self, filename=None):    
        '''Extract basis information from the vasprun.xml'''
        if filename is None:
            if check_exist(self.prefix + ".scf.out"):
                filename = self.prefix + ".scf.out"
            else:
                assert 0, "Cannot find any prefix.scf.out file"
                
        data = qe_io.read_pw_output(filename)
        self.nelec = data['nelec']
        self.nbands = data['nbands']
        self.soc = data['soc']
        self.kpts = data['kpts'][:,:3]
        self.kpts_weight = data['kpts'][:,3]
        self.nkpts = self.kpts.shape[0] 
        self.band = data['eigenvals']
        self.ispin = self.band.shape[0]
        self.atom  = data['atom']
        self.natom  = len(self.atom)
        self.element = [atm[0] for atm in data['species']]
        
        # Make a cell object in the spglib format
        alat = data['alat'] * const.AUTOA
        lattice = data['crystal axes'] * alat   
        positions = data['atom_position']
        numbers = cell_utils.convert_atomtype(self.atom)
        self.cell_init = (lattice, positions, numbers)
        self.cell = self.cell_init

############ Plotting #################

    def get_efermi(self, data=None):
        '''Get Fermi energy'''
        if data is None:
            if check_exist(self.prefix + ".band.out"):
                filename = self.prefix + ".band.out"
            elif check_exist(self.prefix + ".scf.out"):
                filename = self.prefix + ".scf.out"
            data = qe_io.read_pw_output(filename)
            
        if data['hoco'] is not None:
            efermi = data['hoco']
        elif data['efermi'] is not None:
            efermi = data['efermi']
        else:
            # Take the HOMO as the E Fermi
            band = data['eigenvals']
            nelec = int(data['nelec']) 
            if band.shape[0] == 1:
                npairs = nelec//2
                efermi = band[0,:,:npairs].max()
            else:  
                nspin, nkpts, nband = band.shape
                band = np.sort(band.flatten())
                nelecs = nelec * nkpts
                efermi = band[:nelecs].max()
        
        return efermi

    def get_band(self, filename=None):
        '''
        if filename is not specified, band structure is obtained by searching for:
            (1) self.prefix + ".band.out"
            (2) self.prefix + ".scf.out"
        '''
        if filename is None:
            if check_exist(self.prefix + ".band.out"):
                filename = self.prefix + ".band.out"
            elif check_exist(self.prefix + ".scf.out"):
                filename = self.prefix + ".scf.out"
            else:
                assert 0, "Cannot find any band structure file"
                
        # Collecting data from pw_output
        data = qe_io.read_pw_output(filename)
        band = data['eigenvals']
        kpts = data['kpts'][:,:3]       # the last column is the kpt weight
        alat = data['alat'] * const.AUTOA       # in angstrom
        lattice = data['crystal axes'] * alat
        recip_lattice = data['reciprocal axes'] * 2*np.pi / alat
        
        # Find absolute kpts
        abs_kpts = kpts.dot(recip_lattice)                  # From fractional to absolute in A^-1 unit
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())

        efermi = self.get_efermi(data)
        
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
        
    def _generate_pband(self, filename=None, spin=0, gradient=False, lm='spd'):
        '''Processing/collecting the projected band data before the plotting function
            
            Note: In QE, each atom has its own set of atomic orbitals to project onto.
                  For example, Fe has 1s, 2s, p, and d. O only has s and p
                  Unlike QE, VASP decomposes the band into a list of atomic orbitals.  
                  For example, both Fe and O have contributations from s, p, d, .... Hence,
                  the projected wf is sparser in VASP than in QE

            proj_wf = [kpts, band, # of orbitals]
                  
            Examples for lm:
                lm = 'Ni:s ; p ; d'             :   three groups: (1) s of Ni ; (2) all p orbitals ; (3) all d orbitals
                lm = ['Ni:s,pd', 'O1:p;O2']     :   two groups: (1) s,p,d of Ni ; (2) all p orbitals of the 1st O  and all otbitals of O2 
                lm = ['Ni1;O', 'N']             :   two groups: (1) the 1st Ni and all the O atoms ; (2) All N atom
 
            if gradient == True: user has to provide a TWO groups of orbitals  
                for example, lm = 'Ni:s ; p' or ['Ni:s,pd', 'O1:p;O2']  
            
        '''      
       
        if filename is None:
            if check_exist(self.prefix + ".projwfc.out"):
                filename = self.prefix + ".projwfc.out"
            else:
                assert 0, "Cannot find any band structure file"
    
        data = qe_io.read_projwfc_output(filename)
        species = data['species']
        l_list = data['l']
        m_list = np.int64(data['m'])
        proj_wf = data['projwfc'][spin] 

        # Create the possible lm list
        lm_data = {'0': ['s'], '1':['pz', 'px', 'py'], '2':['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']}
        lm_list = []
        for i, l in enumerate(l_list):
            lm_list.append(lm_data[l][m_list[i] - 1])
        
        # Generate pband
        formatted_atom, formatted_lm = str_format.general_lm(lm)
        
        if gradient:        
            assert len(formatted_atom) == 2, "For the gradient plot, you only need to provide two groups of orbitals, for example, lm = 's ; Ni:d' or lm = ['Ni:s', 'O']"

        # Calculate total band and remove zero values
        total = proj_wf.sum(axis=2)
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
        
    def _generate_dos(self, prefix=None, efermi=None, spin=0, lm=None):
        '''Processing/collecting the DOS data before the plotting function
            
            TDOS dimensions: [spin , [E(eV), tdos(E)]]
            
            spin            : spin of DOS.
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        if lm is None: 
            lm = [atom+':s,p,d' for atom in self.element] 
            
        if prefix is None: prefix = self.prefix
        tdos_file = prefix + ".dos" 
        assert check_exist(tdos_file), "Cannot find " + tdos_file
         
        # Compute pDOS
        if check_exist(prefix + ".pdos_tot"):
            # Get total DOS
            total_pdos_data = qe_io.read_pdos_output(prefix + ".pdos_tot")
            tdos = total_pdos_data[spin,:,:2]
            if efermi is None: 
                efermi = self.get_efermi()
            
            # Collect the pdos files
            data = qe_io.read_projwfc_output(prefix + ".projwfc.out")
            species = data['species']
            wfc_id = np.int64(data['wfc'])
            l_list = data['l']
            m_list = np.int64(data['m'])
 
            pdos_data = []
            for i, atm in enumerate(self.atom):
                idx_atom = [j for j, atom in enumerate(species) if atom == atm]
                wfc_idx = np.unique(wfc_id[idx_atom]) 
                for wfc in wfc_idx:
                    filename = prefix + ".pdos_atm#" + str(i + 1) + "(" + atm + ")" + "_wfc#" + str(wfc)
                    if check_exist(filename + "(s)"): 
                        filename = filename + "(s)"
                    elif check_exist(filename + "(p)"):
                        filename = filename + "(p)"                     
                    elif check_exist(filename + "(d)"):
                        filename = filename + "(d)"
                        
                    lm_pdos_data = qe_io.read_pdos_output(filename)[spin] 
                    pdos_data.append(lm_pdos_data[:,2:])
                    
            pdos_data = np.concatenate(pdos_data, axis=1)
            
            # Create the possible lm list
            lm_data = {'0': ['s'], '1':['pz', 'px', 'py'], '2':['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']}
            lm_list = []
            for i, l in enumerate(l_list):
                lm_list.append(lm_data[l][m_list[i] - 1])
                
            formatted_atom, formatted_lm = str_format.general_lm(lm)
            # Compute pDOS
            pdos = []             
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
                    
                    proj_val += (pdos_data[:,idx_lm]).sum(axis=1)
                pdos.append(proj_val)
            pdos = np.asarray(pdos).T 
        else:
            # Get total DOS
            tdos_data = qe_io.read_tdos_output(tdos_file)      
            tdos = tdos_data['dos'][spin,:,:2]
            if efermi is None: 
                efermi = tdos_data['efermi']
            else:
                efermi = 0
            pdos = None 
            
            
        # Shift the energy 
        tdos[:,0] = tdos[:,0] - efermi 
          
        return tdos, pdos
        
    def _generate_kdos(self, prefix=None, efermi=None, spin=0, lm=None, klabel=None):
        '''Processing/collecting the k-resolved DOS data before the plotting function
           The kDOS will be summed over all the lm  
            
            kDOS dimensions: [spin , kpts, [E(eV), tdos(E)]]
            
            spin            : spin of DOS.
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        
        if prefix is None: prefix = self.prefix
        assert check_exist(prefix + ".pdos_tot"), "Cannot find " + tdos_file
        
        formatted_atom, formatted_lm = str_format.general_lm(lm)
        assert len(formatted_atom) == 1, "For kDOS plot, you only need to provide one groups of orbitals, for example, lm = 'Ni:sp', meaning "
            
        # Get total kDOS
        total_kdos_data = qe_io.read_kdos_output(prefix + ".pdos_tot")
        tdos = total_kdos_data[spin,:,:,:2]
        if efermi is None: 
            efermi = self.get_efermi()
            
        # Collect the pdos files
        projwfc_data = qe_io.read_projwfc_output(prefix + ".projwfc.out")
        site = projwfc_data['site']
        species = projwfc_data['species']
        wfc_id = np.int64(projwfc_data['wfc'])
        l_list = projwfc_data['l']
        m_list = np.int64(projwfc_data['m'])
     
        pdos_data = []
        for i, atm in enumerate(self.atom):
            idx_atom = [j for j, atom in enumerate(species) if atom == atm]
            wfc_idx = np.unique(wfc_id[idx_atom]) 
            for wfc in wfc_idx:
                filename = prefix + ".pdos_atm#" + str(i + 1) + "(" + atm + ")" + "_wfc#" + str(wfc)
                if check_exist(filename + "(s)"): 
                    filename = filename + "(s)"
                elif check_exist(filename + "(p)"):
                    filename = filename + "(p)"                     
                elif check_exist(filename + "(d)"):
                    filename = filename + "(d)"
                    
                lm_pdos_data = qe_io.read_kdos_output(filename)[spin] 
                pdos_data.append(lm_pdos_data[:,:,2:])
                
        pdos_data = np.concatenate(pdos_data, axis=2)
        
        # Create the possible lm list
        lm_data = {'0': ['s'], '1':['pz', 'px', 'py'], '2':['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']}
        lm_list = []
        for i, l in enumerate(l_list):
            lm_list.append(lm_data[l][m_list[i] - 1])
                
        # Compute kDOS          
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
                
                proj_val += (pdos_data[:,:,idx_lm]).sum(axis=2)
        pdos = proj_val
            
        # Shift the energy 
        tdos[:,:,0] = tdos[:,:,0] - efermi 
          
        # Compute the kpath projected on 1D
        kpts = projwfc_data['kpts']
        lattice = self.cell[0]
        recip_lattice =  2 * np.pi * np.linalg.inv(lattice).T
        abs_kpts = kpts.dot(recip_lattice)                  # From fractional to absolute in A^-1 unit
        temp_kpts = np.empty_like(abs_kpts)
        temp_kpts[0] = abs_kpts[0]
        temp_kpts[1:] = abs_kpts[:-1] 
        proj_kpath = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
        
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
          
        return tdos, pdos, proj_kpath, sym_kpoint_coor, klabel 