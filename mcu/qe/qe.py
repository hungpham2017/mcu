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
from ..utils import plot
from ..utils.misc import check_exist
from ..vasp import const
from ..cell import utils as cell_utils
from . import qe_io

        
class main:
    def __init__(self,  prefix=None):
        '''
        
        '''
        if prefix is None:
            print("Provide a prefix name for your project, for example, prefix.scf.out, prefix.band.out, etc. Otherwise you would have to specify each output file when performing analysis")
        self.prefix = prefix
        self.get_info() 
        
############ General #################
        
    def get_info(self):    
        '''Extract basis information from the vasprun.xml'''

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
 
    def _generate_band(self, filename=None, efermi=None, spin=0, label=None):
        '''Processing/collecting the band data before the plotting function
           label            : a list of labels and corresponding coordinates
        '''
        band, kpath_frac, proj_kpath, recip_lattice, efermi_ = self.get_band(filename)
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
        # Find absolute coordinates for high symmetric kpoints  
        sym_kpoint_coor = None
        if label is not None:
            assert isinstance(label,list)           # label needs to be a list of labels and corresponding coordinates
            temp = []
            coor_kpts = [] 
            for kpt in label:
                temp.append(kpt[0])
                coor_kpts.append(kpt[1:])
            label = temp       
            coor_kpts = np.asarray(coor_kpts)
            abs_kpts = coor_kpts.dot(recip_lattice)   
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
                
        return band, proj_kpath, sym_kpoint_coor, label, True # This True is just for compariable with VASP
        
    def plot_band(self, efermi=None, label=None, spin=0, save=False, band_color=['#007acc','#808080','#808080'],
                    figsize=(6,6), figname='BAND', xlim=None, ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                                  0 or 1 for spin polarized
                color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                                  
                                  
        '''
        assert isinstance(band_color,list)
        assert len(band_color) == 3
        plot.plot_band(self, efermi=efermi, spin=spin, save=save, band_color=band_color,
                figsize=figsize, figname=figname, xlim=xlim, ylim=ylim, fontsize=fontsize, dpi=dpi, format=format, label=label)
                
    def _generate_pband(self, filename=None, spin=0, style=1, lm='spd'):
        '''Processing/collecting the projected band data before the plotting function
            
            Note: In QE, each atom has its own set of atomic orbitals to project onto.
                  For example, Fe has 1s, 2s, p, and d. O only has s and p
                  Unlike QE, VASP decomposes the band into a common list of atomic orbitals.  
                  For example, both Fe and O have contributations from s, p, d, .... Hence,
                  the projected wf is sparser in VASP than in QE

                  proj_wf = [kpts, band, # of orbitals]
            
            style = 1   : all atoms are considered
                         lm = 's', 'py', 'pz', 'px', 'dxy', 'dyz','dz2','dxz','dx2-y2' or a list of them
                             'sp', 'pd', 'sd', 'spd'  => shortcut
                             each color is used for each lm
                             the marker's radius is proportional to the % of lm 
            style = 2   : considering only a list of orbitals
                         e.g. orb = ['Ni:s','C:pz']
            style = 3   : gradient map to show the character transition
                         lm = 'sp', 'pd', 'sd'                
        '''      
       
        if filename is None:
            if check_exist(self.prefix + ".projwfc.out"):
                filename = self.prefix + ".projwfc.out"
            else:
                assert 0, "Cannot find any band structure file"
    
        data = qe_io.read_projwfc_output(filename)
        site = data['site']
        species = data['species']
        wfc_id = data['wfc']
        l_list = data['l']
        m_list = np.int64(data['m'])
        proj_wf = data['projwfc'][spin] 

        # Create the possible lm list
        lm_data = {'0': ['s'], '1':['pz', 'px', 'py'], '2':['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']}
        lm_list = []
        for i, l in enumerate(l_list):
            lm_list.append(lm_data[l][m_list[i] - 1])
        
        irred_lm_list = list(dict.fromkeys(lm_list))
        
        if style == 1:
            lm_shortcut = ['p','d','sp','ps','pd','dp','sd','ds','spd','sdp','psd','pds','dsp','dps']
            # Check if the lm value is appropriate
            if isinstance(lm,str):
                if lm not in irred_lm_list and lm not in lm_shortcut:
                    raise Exception("WARNING:", lm, "is not recognizable. lm must be", irred_lm_list, lm_shortcut)
                else:
                    if lm == 'p': 
                        lm = [['px','py','pz']]
                    elif lm == 'd': 
                        lm = [['dxy', 'dyz','dz2','dxz','dx2-y2']]
                    elif lm == 'sp': 
                        lm = ['s',['px','py','pz']]
                    elif lm == 'ps': 
                        lm = [['px','py','pz'],'s']
                    elif lm == 'sd': 
                        lm = ['s',['dxy', 'dyz','dz2','dxz','dx2-y2']]
                    elif lm == 'ds': 
                        lm = [['dxy', 'dyz','dz2','dxz','dx2-y2'],'s']
                    elif lm == 'pd': 
                        lm = [['px','py','pz'],['dxy', 'dyz','dz2','dxz','dx2-y2']]
                    elif lm == 'dp': 
                        lm = [['dxy', 'dyz','dz2','dxz','dx2-y2'],['px','py','pz']]
                    elif lm == 'spd':                         
                        lm = ['s',['px','py','pz'],['dxy', 'dyz','dz2','dxz','dx2-y2']]  
                    elif lm == 'sdp':                         
                        lm = ['s',['dxy', 'dyz','dz2','dxz','dx2-y2'],['px','py','pz']] 
                    elif lm == 'psd':                         
                        lm = [['px','py','pz'],'s',['dxy', 'dyz','dz2','dxz','dx2-y2']]   
                    elif lm == 'pds':                         
                        lm = [['px','py','pz'],['dxy', 'dyz','dz2','dxz','dx2-y2'],'s']  
                    elif lm == 'dsp':                         
                        lm = [['dxy', 'dyz','dz2','dxz','dx2-y2'],'s',['px','py','pz']] 
                    elif lm == 'dps':                         
                        lm = [['dxy', 'dyz','dz2','dxz','dx2-y2'],['px','py','pz'],'s']                          
                    else:
                        lm = [lm]
            elif isinstance(lm,list):
                for each_lm in lm:
                    if isinstance(each_lm,str):
                        if each_lm not in irred_lm_list:
                            raise Exception("WARNING:", lm, "is not recognizable. lm must be one of these", irred_lm_list)
                    else:
                        for orb in each_lm:
                            if orb not in irred_lm_list:
                                raise Exception("WARNING:", orb , "is not recognizable. lm must be one of these", irred_lm_list)                        
            else:
                raise Exception("lm is not recognizable")
                
            # Compute pband
            total = proj_wf.sum(axis=2)       # Sum over the orbitals --> [kpt,band]
            shape = total.shape
            idx_zeros = total.flatten() < 0.0001
            total = total.flatten()
            total[idx_zeros] = 1.0
            total = total.reshape(shape)
            
            pband = [] 
            for each_lm in lm:
                if isinstance(each_lm,str):  
                    idx_lm = [i for i, x in enumerate(lm_list) if x in each_lm]
                    proj_val = (proj_wf[:,:,idx_lm]).sum(axis=2) / total
                else:
                    proj_val = 0
                    for orb in each_lm:
                        idx_lm = [i for i, x in enumerate(lm_list) if x in orb]
                        proj_val += (proj_wf[:,:,idx_lm]).sum(axis=2) / total
                pband.append(proj_val)
 
            pband = np.asarray(pband)
          
        elif style == 2:
            if isinstance(lm,str):
                atom, lm_ = lm.split(':')
                lm_  = lm_.split(',') 
                temp = []
                for i in range(len(lm_)):               
                    if lm_[i] == 'p': 
                        for m in ['px','py','pz']: temp.append(m)
                    elif lm_[i] == 'd': 
                        for m in ['dxy', 'dyz','dz2','dxz','dx2-y2']: temp.append(m)
                    else:
                        temp.append(lm_[i])
                lms = [temp]
                atoms = [atom]
            elif isinstance(lm,list):
                atoms = []
                lms = []   
                for orb in lm:
                    atom, lm_ = orb.split(':')
                    lm_  = lm_.split(',') 
                    temp = []
                    for i in range(len(lm_)):               
                        if lm_[i] == 'p': 
                            for m in ['px','py','pz']: temp.append(m)
                        elif lm_[i] == 'd': 
                            for m in ['dxy', 'dyz','dz2','dxz','dx2-y2']: temp.append(m)
                        else:
                            temp.append(lm_[i])
                    atoms.append(atom)
                    lms.append(temp)
  
            # Compute pband
            total = proj_wf.sum(axis=2)       # Sum over the orbitals --> [kpt,band]
            shape = total.shape
            idx_zeros = total.flatten() < 0.0001
            total = total.flatten()
            total[idx_zeros] = 1.0
            total = total.reshape(shape)
            
            pband = [] 
            for i, atm in enumerate(atoms):
                idx_atom = [j for j, atom in enumerate(species) if atom == atm]
                idx_lm = [idx for idx in idx_atom if lm_list[idx] in lms[i]]
                proj_val = (proj_wf[:,:,idx_lm]).sum(axis=2)
                pband.append(proj_val/total)
            pband = np.asarray(pband)
            
        elif style == 3: 
            lm_shortcut = ['sp', 'sd', 'pd']
            if isinstance(lm,str):
                if lm not in lm_shortcut:
                    raise Exception("WARNING:", lm, "is not recognizable. lm must be", lm_shortcut)
                else:
                    if lm == 'sp': 
                        lm = ['s',['px','py','pz']]
                    elif lm == 'sd': 
                        lm = ['s',['dxy', 'dyz','dz2','dxz','dx2-y2']]
                    elif lm == 'pd': 
                        lm = [['px','py','pz'],['dxy', 'dyz','dz2','dxz','dx2-y2']]
                    else:
                        raise Exception("WARNING:", lm, "is not recognizable. lm must be one of these", lm_shortcut, "or a list")
            elif isinstance(lm,list):
                assert len(lm) == 2          # Only two orbital 
                for each_lm in lm:
                    if isinstance(each_lm,str):
                        if each_lm not in irred_lm_list:
                            raise Exception("WARNING:", lm, "is not recognizable. lm must be one of these", irred_lm_list)
                    else:
                        for orb in each_lm:
                            if orb not in irred_lm_list:
                                raise Exception("WARNING:", orb , "is not recognizable. lm must be one of these", irred_lm_list) 
            else:
                raise Exception("lm is not recognizable")
                
            # Compute pband
            pband = [] 
            for each_lm in lm:                  # only two lm
                if isinstance(each_lm,str):  
                    idx_lm = [i for i, x in enumerate(lm_list) if x in each_lm]
                    proj_val = (proj_wf[:,:,idx_lm]).sum(axis=2)
                else:
                    proj_val = 0
                    for orb in each_lm:
                        idx_lm = [i for i, x in enumerate(lm_list) if x in orb]
                        proj_val += (proj_wf[:,:,idx_lm]).sum(axis=2)
                pband.append(proj_val)
            pband = np.asarray(pband)
            pband = pband[0]/(pband.sum(axis=0))
        else:
            raise Exception('mcu currently supports only style: 0,1,2')
        
        return pband
                
                    
    def plot_pband(self, efermi=None, spin=0, label=None, style=1, lm='spd', band=None, color=None, band_color=['#007acc','#808080','#808080'],
                    scale=1.0, alpha=0.5, cmap='bwr', edgecolor='none', facecolor=None, marker=None,
                    legend=None, loc="upper right", legend_size=1.0,
                    save=False, figname='pBAND', figsize=(6,6), xlim=None, ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot projected band structure
           Please see mcu.utils.plot.plot_pband for full documents        
        '''
        plot.plot_pband(self, efermi=efermi, spin=spin, label=label, style=style, lm=lm, band=band, color=color, band_color=band_color,
                    scale=scale, alpha=alpha, cmap=cmap, edgecolor=edgecolor, facecolor=facecolor, marker=marker,
                    legend=legend, loc=loc, legend_size=legend_size,
                    save=save, figname=figname, figsize=figsize, xlim=xlim, ylim=ylim, fontsize=fontsize, dpi=dpi, format=format)

    def _generate_dos(self, prefix=None, efermi=None, spin=0, lm=None):
        '''Processing/collecting the DOS data before the plotting function
            
            TDOS dimensions: [spin , [E(eV), tdos(E)]]
            
            spin            : spin of DOS.
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        
        if prefix is None: prefix = self.prefix
        tdos_file = prefix + ".dos" 
        assert check_exist(tdos_file), "Cannot find " + tdos_file
         
        # Compute pDOS
        if check_exist(prefix + ".pdos_tot"):
            # Collecting group of lm
            if isinstance(lm,str):
                atom, lm_ = lm.split(':')
                lm_  = lm_.split(',') 
                temp = []
                for i in range(len(lm_)):               
                    if lm_[i] == 'p': 
                        for m in ['px','py','pz']: temp.append(m)
                    elif lm_[i] == 'd': 
                        for m in ['dxy', 'dyz','dz2','dxz','dx2-y2']: temp.append(m)
                    else:
                        temp.append(lm_[i])
                lms = [temp]
                atoms = [atom]
                
            elif isinstance(lm,list):
                atoms = []
                lms = []   
                for orb in lm:
                    atom, lm_ = orb.split(':')
                    lm_  = lm_.split(',') 
                    temp = []
                    for i in range(len(lm_)):               
                        if lm_[i] == 'p': 
                            for m in ['px','py','pz']: temp.append(m)
                        elif lm_[i] == 'd': 
                            for m in ['dxy', 'dyz','dz2','dxz','dx2-y2']: temp.append(m)
                        else:
                            temp.append(lm_[i])
                    atoms.append(atom)
                    lms.append(temp)
            
            # Get total DOS
            total_pdos_data = qe_io.read_pdos_output(prefix + ".pdos_tot")
            tdos = total_pdos_data[spin,:,:2]
            if efermi is None: 
                efermi = self.get_efermi()
            
            # Collect the pdos files
            data = qe_io.read_projwfc_output(prefix + ".projwfc.out")
            site = data['site']
            species = data['species']
            wfc_id = np.int64(data['wfc'])
            l_list = data['l']
            m_list = np.int64(data['m'])
 
            # Create the possible lm list
            lm_data = {'0': ['s'], '1':['pz', 'px', 'py'], '2':['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']}
            lm_list = []
            for i, l in enumerate(l_list):
                lm_list.append(lm_data[l][m_list[i] - 1])
            
            irred_lm_list = list(dict.fromkeys(lm_list))
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
                    
            pdos_data = np.hstack(pdos_data)
            
            # Compute pDOS
            pdos = [] 
            for i, atm in enumerate(atoms):
                idx_atom = [j for j, atom in enumerate(species) if atom == atm]
                idx_lm = [idx for idx in idx_atom if lm_list[idx] in lms[i]]
                proj_val = (pdos_data[:,idx_lm]).sum(axis=1)
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
        
    def plot_dos(self, style=1, efermi=None, spin=0, lm=None, color=None,
                    legend=None, loc="upper right", fill=True, alpha=0.2,
                    save=False, figname='DOS', figsize=(6,3), elim=(-6,6), yscale=1.1, fontsize=18, dpi=600, format='png'):
        '''Plot projected band structure
           Please see mcu.utils.plot.plot_dos for full documents 
        '''
        plot.plot_dos(self, style=style, efermi=efermi, spin=spin, lm=lm, color=color,
                legend=legend, loc=loc, fill=fill, alpha=alpha,
                save=save, figname=figname, figsize=figsize, elim=elim, yscale=yscale, fontsize=fontsize, dpi=dpi, format=format)
        
    def _generate_kdos(self, prefix=None, efermi=None, spin=0, lm=None):
        '''Processing/collecting the k-resolved DOS data before the plotting function
            
            kDOS dimensions: [spin , kpts, [E(eV), tdos(E)]]
            
            spin            : spin of DOS.
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        
        if prefix is None: prefix = self.prefix
        tdos_file = prefix + ".dos" 
        assert check_exist(tdos_file), "Cannot find " + tdos_file
         
        # Compute pDOS
        if check_exist(prefix + ".pdos_tot"):
            # Collecting group of lm
            if isinstance(lm,str):
                atom, lm_ = lm.split(':')
                lm_  = lm_.split(',') 
                temp = []
                for i in range(len(lm_)):               
                    if lm_[i] == 'p': 
                        for m in ['px','py','pz']: temp.append(m)
                    elif lm_[i] == 'd': 
                        for m in ['dxy', 'dyz','dz2','dxz','dx2-y2']: temp.append(m)
                    else:
                        temp.append(lm_[i])
                lms = [temp]
                atoms = [atom]
                
            elif isinstance(lm,list):
                atoms = []
                lms = []   
                for orb in lm:
                    atom, lm_ = orb.split(':')
                    lm_  = lm_.split(',') 
                    temp = []
                    for i in range(len(lm_)):               
                        if lm_[i] == 'p': 
                            for m in ['px','py','pz']: temp.append(m)
                        elif lm_[i] == 'd': 
                            for m in ['dxy', 'dyz','dz2','dxz','dx2-y2']: temp.append(m)
                        else:
                            temp.append(lm_[i])
                    atoms.append(atom)
                    lms.append(temp)
            
            # Get total DOS
            total_pdos_data = qe_io.read_pdos_output(prefix + ".pdos_tot")
            tdos = total_pdos_data[spin,:,:2]
            if efermi is None: 
                efermi = self.get_efermi()
            
            # Collect the pdos files
            data = qe_io.read_projwfc_output(prefix + ".projwfc.out")
            site = data['site']
            species = data['species']
            wfc_id = np.int64(data['wfc'])
            l_list = data['l']
            m_list = np.int64(data['m'])
 
            # Create the possible lm list
            lm_data = {'0': ['s'], '1':['pz', 'px', 'py'], '2':['dz2', 'dxz', 'dyz', 'dx2-y2', 'dxy']}
            lm_list = []
            for i, l in enumerate(l_list):
                lm_list.append(lm_data[l][m_list[i] - 1])
            
            irred_lm_list = list(dict.fromkeys(lm_list))
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
                    
            pdos_data = np.hstack(pdos_data)
            
            # Compute pDOS
            pdos = [] 
            for i, atm in enumerate(atoms):
                idx_atom = [j for j, atom in enumerate(species) if atom == atm]
                idx_lm = [idx for idx in idx_atom if lm_list[idx] in lms[i]]
                proj_val = (pdos_data[:,idx_lm]).sum(axis=1)
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