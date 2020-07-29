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
from . import crystal_io
from . import utils as crystal_utils

        
class main(cell.main, plot.main):
    def __init__(self,  prefix="outfile"):
        assert prefix is not None, "Provide a prefix name for your project, for example, prefix.scf.out, prefix.band.out, etc."
        self.prefix = prefix
        self.get_info() 
        
############ General #################
    def get_info(self, filename=None):    
        '''Extract basis information from the vasprun.xml'''
        if filename is None:
            if check_exist(self.prefix + ".out"):
                filename = self.prefix + ".out"
            else:
                assert 0, "Cannot find any prefix.out file"
                
        data = crystal_io.read_out(filename)
        self.nelec = data['nelec']
        self.nelec_core = data['nelec core']
        self.nao = data['nao']
        self.basis = data['basis']
        self.atom  = data['atom']
        self.natom  = len(self.atom)
        self.element = list(dict.fromkeys(self.atom ))
        self.kpts = data['kpts']
        
        # Make a cell object in the spglib format
        lattice = data['lattice']
        positions = data['atom_position']
        numbers = cell_utils.convert_atomtype(self.atom)
        self.cell_init = (lattice, positions, numbers)
        self.cell = self.cell_init
        
############ Preprocessing ################# 
    def make_DOS_input(self, lm='None', filename=None):
        '''Giving a lm string, export a d3 file for DOS calculation.
        '''
        if lm is None: lm = 'spdf'
        if filename is None: filename = self.prefix + '.d3'
        formatted_atom, formatted_lm = str_format.general_lm(lm)
        
        pDOS_orb_idx = []             
        for i, atoms in enumerate(formatted_atom):  
            proj_orbs = []
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
                            elif count // nwfc == id - 1: 
                                idx_atom.append(n)
                            count += 1 
 
                if atom is None: 
                    atom_ = ""
                    available_orbs = crystal_utils.basis_short2long(self.basis)
                    available_orbs = sum(list(available_orbs.values()),[])
                else:
                    atom_, id = str_format.format_atom(atom)
                    available_orbs = crystal_utils.basis_short2long(self.basis)[atom_]
                
                available_orbs = list(dict.fromkeys(['None'] + available_orbs))
                
                for atm_idx in idx_atom:
                    for orb in formatted_lm[i][j]:
                        assert str(orb) in available_orbs, "This is wrong: " + str(orb) + ". Check the lm string/list. Available basis functions " + atom_ +  " are: " +  ", ".join(available_orbs)
                        idx = crystal_utils.orb_index(self.atom, self.basis, atm_idx, orb)
                        proj_orbs.append(idx)
            proj_orbs = sum(proj_orbs,[])
            proj_orbs.sort()
            pDOS_orb_idx.append(proj_orbs)

        # Export d3 file using pDOS_orb_idx list
        nDOS = len(pDOS_orb_idx)
        npoints = 1000
        first_band = self.nelec_core//2 + 1
        last_band = self.nao
        store_dos = 1
        legendre = 14
        printing = 0
        with open(filename, 'w') as f:
            f.write('NEWK\n')
            f.write('12 12\n')
            f.write('1 0\n')
            f.write('DOSS\n')  
            f.write('%d %d %d %d %d %d %d\n' % (nDOS, npoints, first_band, last_band, store_dos, legendre, printing))   
            for dos in pDOS_orb_idx:
                norb = len(dos)
                orb_string = [str(orb) for orb in dos]
                f.write('%d %s\n' % (norb, ' '.join(orb_string)))         
            f.write('END\n') 

        
############ Plotting ################# 
    def get_band(self, filename=None, phonon=False, gamma_correct=False, threshold=8.06554):
        '''Make a band from from 
           NOTE: the proj_kpath is computed from the dk and nkp. Due to the round out error in f25, the computed high symmetric k-point coordinates won't be exactly the same as values obtained from *.BAND file.
           
           if phonon == True: the phonon band is read. Three small negative acoustic modes at Gamma can be removed by using the gamma_correct keyword and threshold 1 meV = 8.06554 cm^-1  
        '''
        if filename is None:
            if check_exist(self.prefix + ".f25"):
                filename = self.prefix + ".f25"
            elif check_exist("fort.25"):
                filename = "fort.25"
                print('Found fort.25. Band structure is extracted from fort.25')            
            else:
                assert 0, "Cannot find " + self.prefix + ".f25" + " or fort.25. Check if you has band structure file"

        data = crystal_io.read_f25(filename)[0]
        assert data is not None, "Cannot find BAND information, check if you have generated fort.25: " + filename
        temp = []
        ihferm_list = []
        proj_kpath = []
        sym_kpoint_coor = [0.0]
        shift = 0
        for block in data:
            ihferm = block['ihferm']
            nband = block['nband']
            nkp = block['nkp']
            dum = block['dum']
            dk = block['dk']
            efermi = block['efermi']
            emin = block['emin']
            emax = block['emax']
            ivalues = block['values']
            eigenvals = block['eigenvals']

            # Compute hight symmetric k-pts coordinates
            temp.append(eigenvals.reshape(-1, nband))
            path = np.arange(nkp)*dk + shift
            proj_kpath.append(path)
            shift = path[-1] + dk
            sym_kpoint_coor.append(path[-1])
        
        if ihferm == 0 or ihferm == 2:
            band = np.float64([np.vstack(temp)])
            proj_kpath = np.hstack(proj_kpath)
        elif ihferm == 1 or ihferm == 3:
            nblock = len(temp) // 2
            band_up = np.vstack(temp[:nblock])            
            band_down = np.vstack(temp[nblock:])  
            band = np.float64([band_up, band_down])
            proj_kpath = np.hstack(proj_kpath[:nblock])
            
        sym_kpoint_coor = np.float64(sym_kpoint_coor)
        
        if phonon:
            if gamma_correct:
                # Should be checked first before using the correction to make sure it is really acoustic phonon modes
                nspin, nkpts, nband = band.shape
                band = band.flatten()
                imag_mode_idx = band < 0.0
                imag_mode = band[imag_mode_idx]
                imag_mode_idx2 = imag_mode > -threshold
                imag_mode[imag_mode_idx2] = 0.0
                band[imag_mode_idx] = imag_mode
                band = band.reshape(nspin, nkpts, nband)
        else:
            band = const.AUTOEV * band
            efermi = const.AUTOEV * efermi

        return band, proj_kpath, sym_kpoint_coor, efermi

    def get_bandgap(self, filename=None, efermi=None):
        '''Get the bandgap'''
        if filename is None:
            if check_exist(self.prefix + ".f25"):
                filename = self.prefix + ".f25"
            elif check_exist("fort.25"):
                filename = "fort.25"
                print('Found fort.25. Band structure is extracted from fort.25')            
            else:
                assert 0, "Cannot find " + self.prefix + ".f25" + " or fort.25. Check if you has band structure file"
                
        band, proj_kpath, sym_kpoint_coor, efermi_ = self.get_band(filename)
        if efermi is None: efermi = efermi_
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
                
                # TODO: the kpath_frac currently cannot be obtained form f25.
                # Other outputs are needed
                # kpath_frac = self.cp2k_io.kpath_frac[set_block]
                # print('  E(VBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (VBM[vbm_idx], 
                                                                # kpath_frac[vbm_idx,0], kpath_frac[vbm_idx,1], kpath_frac[vbm_idx,2]))
                # print('  E(CBM) = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (CBM[cbm_idx], 
                                                                # kpath_frac[cbm_idx,0], kpath_frac[cbm_idx,1], kpath_frac[cbm_idx,2]))
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
        '''
        if filename is None:
            if check_exist(self.prefix + ".f25"):
                filename = self.prefix + ".f25"
            elif check_exist("fort.25"):
                filename = "fort.25"
                print('Found fort.25. Band structure is extracted from fort.25')            
            else:
                assert 0, "Cannot find " + self.prefix + ".f25" + " or fort.25. Check if you has band structure file"
                
        band, proj_kpath, sym_kpoint_coor, efermi_ = self.get_band(filename)
        if efermi is None: efermi = efermi_  
        band = band[spin] - efermi
        
        # Process the labels for k-point
        if klabel is not None:
            klabel, coor_kpts = str_format.format_klabel(klabel) 
            assert len(klabel) == len(sym_kpoint_coor), "The number of k label must be " + str(len(sym_kpoint_coor))

        return band, proj_kpath, sym_kpoint_coor, klabel
        
    def _generate_phononband(self, unit="CM", gamma_correct=False, threshold=8.06554, spin=0, klabel=None):
        '''Processing/collecting the phnon band (cm^-1) data before the plotting function
           Unit: CM = CM^-1, THZ = THz, MEV = meV
        '''
        band, proj_kpath, sym_kpoint_coor, efermi_ = self.get_band(phonon=True, gamma_correct=gamma_correct, threshold=threshold)        
        if unit.lower() == "thz":
            band = const.CMTOTHZ * band
        elif unit.lower() == "mev":
            band = const.CMTOMEV * band
            
        # Process the labels for k-point
        if klabel is not None:
            klabel, coor_kpts = str_format.format_klabel(klabel) 
            assert len(klabel) == len(sym_kpoint_coor), "The number of k label must be " + str(len(sym_kpoint_coor))
            
        return band[spin], proj_kpath, sym_kpoint_coor, klabel
        
    def _generate_dos(self, filename=None, efermi=None, spin=0, lm=None):
        '''Processing/collecting the DOS data before the plotting function
            
            TDOS dimensions: [spin , [E(eV), tdos(E)]]
            
            spin            : spin of DOS.
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        if lm is None: 
            lm = [atom+':s,p,d' for atom in self.element] 
        
        if filename is None:
            if check_exist(self.prefix + ".f25"):
                filename = self.prefix + ".f25"
            elif check_exist("fort.25"):
                filename = "fort.25"
                print('Found fort.25. DOS is extracted from fort.25')            
            else:
                assert 0, "Cannot find " + self.prefix + ".f25" + " or fort.25. Check if you has DOS file"

        # Compute pDOS
        data = crystal_io.read_f25(filename)[1]
        assert data is not None, "Cannot find DOS information, check if you have generated fort.25: " + filename

        dos_data = []
        for block in data:
            ihferm = block['ihferm']
            nrow = block['nrow']
            ncol = block['ncol']
            dx = block['dx']
            dy = block['dy']
            efermi = const.AUTOEV * block['efermi']
            emin = block['emin']
            emax = block['emax']
            ivalues = block['values']
            dos_data.append(block['dos'])

        if ihferm == 0 or ihferm == 2:
            band = np.float64([np.vstack(temp)])
            proj_kpath = np.hstack(proj_kpath)
        elif ihferm == 1 or ihferm == 3:
            nblock = len(temp) // 2
            band_up = np.vstack(temp[:nblock])            
            band_down = np.vstack(temp[nblock:])  
            band = np.float64([band_up, band_down])
            proj_kpath = np.hstack(proj_kpath[:nblock])

        # Compute energy points:
        energy_path = const.AUTOEV * (np.arange(ncol)*dy + emax) - efermi
        tdos = np.float64([energy_path, dos_data[-1]]).T
        pdos = None

        if len(dos_data) != 1:
            pdos = np.asarray(dos_data[:-1]).T
          
        return tdos, pdos       