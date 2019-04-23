#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import mcu
from mcu.vasp import utils, read
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
        
class VASP:
    def __init__(self, path='./', vaspruns='vasprun', outcars='OUTCAR'):
        '''
            path        : the project directory
            vaspruns    : a str or a list of string as names for *.xml files
            outcars     : a str or a list of string as names for OUTCAR files
        '''
        
        # Create vasprun object(s)
        if isinstance(vaspruns, str):                   # For one vasprun.xml file    
            self.vasprun = read.vasprun(path + '/' + vaspruns + '.xml')
            self.useOUTCAR = False
            self.outcar = read.OUTCAR(path + '/' + outcars)
            if self.outcar.success == True: self.useOUTCAR = True
                
            self.get_info(self.vasprun)
        elif isinstance(vaspruns, list):                # For multiple vasprun.xml file
            self.vasprun = []
            for xml in vaspruns:
                xml_file = path + '/' + xml + '.xml'
                if not utils.check_exist(xml_file):
                    print('Cannot find:', xml_file)
                    break
                self.vasprun.append(read.vasprun(xml_file))

            self.useOUTCAR = False
            if isinstance(outcars, list):               
                assert len(outcars) == len(vaspruns)
                self.outcar = []
                for outcar in outcars:
                    outcar_file = path + '/' + outcar
                    if not utils.check_exist(outcar_file):
                        print('Cannot find:', outcar_file)
                        break
                    self.outcar.append(read.OUTCAR(outcar_file))  
                    self.useOUTCAR = True                    
            self.get_info(self.vasprun[0])      # Only get info for the first vasprun.xml
        else:
            print('Provide a string or a list of names for *.xml file')
    
    def get_info(self, vasprun):    
        '''Extract basis information from the vasprun.xml'''
        
        electronic = vasprun.parameters['electronic']
        self.nelec = electronic.general['NELECT']    
        self.nbands = electronic.general['NBANDS']
        self.lsorbit = electronic.spin['LSORBIT']
        self.ispin = electronic.spin['ISPIN']    
        self.kpts = vasprun.kpoints['kpointlist']
        self.nkpts = self.kpts.shape[0]        
        self.natom  = vasprun.natom 
        self.atom  = vasprun.atom     
        self.atm  = vasprun.atm 
        self.atype = [atom[1] for atom in vasprun.types]
        self.get_efermi()
        
    def get_efermi(self):
        '''Extract E_fermi either from vasprun.xml or OUTCAR'''
        if isinstance(self.vasprun, mcu.vasp.read.vasprun):
            self.vasprun.get_dos()
            if hasattr(self.vasprun,'efermi'):
                self.efermi = self.vasprun.efermi
            else:
                if self.useOUTCAR == False:
                    print ("Fermi level need to be read from OUTCAR")
                else:
                    self.efermi = self.outcar.efermi 
        elif isinstance(self.vasprun, list):        
            self.efermi = []
            for i in range(len(self.vasprun)):
                self.vasprun[i].get_dos()
                if hasattr(self.vasprun[i],'efermi'):
                    self.efermi.append(self.vasprun[i].efermi)
                else:
                    if self.useOUTCAR == False:
                        print ("Fermi level need to be read from OUTCAR")
                    else:
                        self.efermi.append(self.outcar[i].efermi)                 
            
    def get_bandgap(self, efermi=None):
        '''Get the bandgap'''
        
        # Get the fermi level
        if efermi == None: efermi = self.efermi
            
        if isinstance(self.vasprun,mcu.vasp.read.vasprun):              # For one vasprun.xml file
            assert isinstance(efermi,float) 
            self.vasprun.get_band()
            self.band = self.vasprun.band[:,:,:,0]
            self.co_occ = self.vasprun.band[:,:,:,1] 
            self.co_occ_ = self.co_occ > 0.5       
            electronic = self.vasprun.parameters['electronic']
            self.kpts = self.vasprun.kpoints['kpointlist']
            nkpts = self.kpts.shape[0]
        elif isinstance(self.vasprun,list):                             # For multiple vasprun.xml file
            assert isinstance(efermi,list)
            for i in range(len(self.vasprun)): 
                assert isinstance(efermi[i],float)
                
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
                for i, vasprun in enumerate(self.vasprun):
                    vasprun.get_band()
                    band = vasprun.band[spin,:,:,0]
                    kpts = vasprun.kpoints['kpointlist']
                    weight = vasprun.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    kpts, band = kpts[nonzero:], band[nonzero:]
                    co_occ = vasprun.band[spin,nonzero:,:,1]
                    co_occ_ = band < efermi[i] 
                    bands = np.vstack([bands,band])
                    kptss = np.vstack([kptss,kpts])
                    co_occ1 = np.vstack([co_occ1,co_occ])
                    co_occ2 = np.vstack([co_occ2,co_occ_])
                band_spin.append(bands[1:])  
                co_occ_spin1.append(co_occ1[1:])   
                co_occ_spin2.append(co_occ2[1:])  
            self.kpts, self.band = np.asarray(kptss[1:]), np.asarray(band_spin)
            nkpts = self.kpts.shape[0]
            self.co_occ, self.co_occ_ = np.asarray(co_occ_spin1), np.asarray(co_occ_spin2)
            
        bandedge = np.zeros([self.ispin,nkpts,2,2])
        self.bandgap = []
        for spin in range(self.ispin):
            print('Spin:', spin)        
            for kpt in range(nkpts):
                band_kpt = self.band[spin,kpt]
                occ = self.co_occ_[spin,kpt]               
                homo_idx = np.count_nonzero(occ) - 1
                lumo_idx = homo_idx + 1               
                bandedge[spin,kpt,0,0] = band_kpt[homo_idx]
                bandedge[spin,kpt,0,1] = self.co_occ[spin,kpt,homo_idx]
                bandedge[spin,kpt,1,0] = band_kpt[lumo_idx]
                bandedge[spin,kpt,1,1] = self.co_occ[spin,kpt,lumo_idx]
                
            vbm_idx = np.argmax(bandedge[spin,:,0,0])
            cbm_idx = np.argmin(bandedge[spin,:,1,0])
            direct = False
            if vbm_idx == cbm_idx: direct = True
            print('  E(VBM) = %7.4f with occ = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (bandedge[spin,vbm_idx,0,0], bandedge[spin,vbm_idx,0,1], 
                                                            self.kpts[vbm_idx,0],self.kpts[vbm_idx,1],self.kpts[vbm_idx,2]))
            print('  E(CBM) = %7.4f with occ = %7.4f at k = [%6.4f,%6.4f,%6.4f]' % (bandedge[spin,cbm_idx,1,0], bandedge[spin,cbm_idx,1,1], 
                                                            self.kpts[cbm_idx,0],self.kpts[cbm_idx,1],self.kpts[cbm_idx,2]))
            bandgap = bandedge[spin,cbm_idx,1,0] - bandedge[spin,vbm_idx,0,0] 
            self.bandgap.append(bandgap)  
            if direct == True: 
                print('  Direct bandgap   : %6.3f' % (bandgap))             
            else:  
                print('  Indirect bandgap : %6.3f' % (bandgap))              
                gap1 = bandedge[spin,cbm_idx,1,0] - bandedge[spin,cbm_idx,0,0]
                gap2 = bandedge[spin,vbm_idx,1,0] - bandedge[spin,vbm_idx,0,0]  
                direct_gap = gap1
                if gap1 > gap2: direct_gap = gap2
                print('  Direct bandgap   : %7.4f' % (direct_gap))                   
                
                
    def _generate_band(self, vasprun, efermi=None, spin=0, label=None):
        '''Processing/collecting the band data before the plotting function
        '''
        
        # Get the fermi level
        if efermi == None: efermi = self.efermi        
        
        sym_kpoint_coor = None
        band = None
        path = None
        conventional = False
            
        if isinstance(vasprun,mcu.vasp.read.vasprun) and vasprun.kpoints['type'] == 1: # For conventional band structure calculation 
            if label != None:
                assert isinstance(label,str)     # label needs to be a string in the format,e.g. 'A-B-C-D'
                label = label.split('-')
            assert isinstance(efermi,float)
            conventional = True  
            vasprun.get_band()
            band = vasprun.band[spin][:,:,0]
            kpts = vasprun.kpoints['kpointlist']
            kpts, band = utils.rm_redundant_band(kpts, band) 
            
            # Find absolute kpts and shift the band
            b = vasprun.cell_final[1]               # Get the reciprocal lattice
            abs_kpts = kpts.dot(b)                  # From fractional to absolute
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            path = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            band = band - efermi               # Efermi is set at 0 eV
            
            highsym_kpt = vasprun.kpoints['points']
            nkpts = highsym_kpt.shape[0]
            sym_kpoint_coor = [0.0]
            for kpt in range(nkpts-2):
                idx = ((path.shape[1] + nkpts - 2)//(nkpts-1) - 1) * (kpt+1)
                coor = path[0,idx]         
                sym_kpoint_coor.append(coor)
            sym_kpoint_coor.append(1.0*path.max())   
            sym_kpoint_coor = np.asarray(sym_kpoint_coor)
                     
        else:
            if isinstance(vasprun,mcu.vasp.read.vasprun):                       # For one vasprun.xml file
                assert isinstance(efermi,float)
                vasprun.get_band()
                band = vasprun.band[spin][:,:,0]
                kpts = vasprun.kpoints['kpointlist']
                if vasprun.kpoints['type'] == 0:
                    weight = vasprun.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    kpts, band = kpts[nonzero:], band[nonzero:]
                band = band - efermi
            elif isinstance(vasprun,list):                                      # For multiple vasprun.xml file
                assert isinstance(efermi,list)
                for i in range(len(self.vasprun)): 
                    assert isinstance(efermi[i],float)
                
                electronic = vasprun[0].parameters['electronic']
                nbands = electronic.general['NBANDS']
                bands = np.zeros([1,nbands])
                kptss = np.zeros([1,3])
                for i, run in enumerate(vasprun):
                    run.get_band()
                    band = run.band[spin][:,:,0]
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
            path = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())

            # Find absolute coordinates for high symmetric kpoints  
            if label != None:
                assert isinstance(label,list)           # label needs to be a list of labels and corresponding coordinates
                temp = []
                coor_kpts = [] 
                for kpt in label:
                    temp.append(kpt[0])
                    coor_kpts.append(kpt[1:])
                label = temp       
                coor_kpts = np.asarray(coor_kpts)
                abs_kpts = coor_kpts.dot(b)   
                temp_kpts = np.empty_like(abs_kpts)
                temp_kpts[0] = abs_kpts[0]
                temp_kpts[1:] = abs_kpts[:-1] 
                sym_kpoint_coor = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 
        
        return band, path, sym_kpoint_coor, label, conventional
                
    def plot_band(self, efermi=None, spin=0, label=None, save=False, band_color=['#007acc','#808080','#808080'],
                    figsize=(6,6), figname='BAND', xlim=None, ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : 0  for spin unpolarized and LSORBIT = .TRUE.
                                  0 or 1 for spin polarized
                label           : label for high symmetric points, e.g. 'G-X-L-W-G'
                                  if hybridXC=True, the lavel should be a list of labels plus their coordinates
                color           : a list of three color codes for band curves, high symmetric kpoint grid, and Fermi level
                                  
                                  
        '''
        
        assert isinstance(band_color,list)
        assert len(band_color) == 3
        
        band, path, sym_kpoint_coor, label, conventional = self._generate_band(self.vasprun, efermi, spin, label)  

        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        yrange = (-50,50)
        
        # Plot the high symmetric kpoint grid
        if conventional == True or label != None:
            for kpt in range(sym_kpoint_coor.shape[0]):
                ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color=band_color[1],linewidth=1.0)

        if label != None and xlim == None:
            nkpts = len(label)
            assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
            for kpt in range(nkpts):   
                point = label[kpt]
                if point == 'G': point = r'$\Gamma$'
                ax.text(sym_kpoint_coor[kpt]/path.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                        color='black', fontsize=fontsize)    

        # Plot bands            
        ax.plot([0,path.max()],[0,0],color=band_color[2],linewidth=1.0, dashes=[6,3])       # Fermi level
        for ith in range(band.shape[1]):
            ax.plot(path.T,band[:,ith],color=band_color[0],linewidth=1.0)    
             
        # Graph adjustments             
        ax.tick_params(labelsize=fontsize)
        if xlim == None:
            plt.xlim([0,path.max()])
            plt.xticks([])
            plt.xlabel('k', size=fontsize+4)
        else:
            plt.xlim(xlim)
            plt.xlabel('k ' + r'($\AA^{-1}$)', size=fontsize+4)
        ax.xaxis.set_label_coords(0.5, -0.08) 
        plt.ylabel('Energy (eV)', size=fontsize+4)        
        plt.ylim(ylim)
        plt.tight_layout()
        if save == True: 
            fig.savefig('Band.'+format,dpi=dpi,format=format)      
        else:
            plt.show()            
            
        
    def _generate_pband(self, vasprun, spin=0, style=1, lm='spd', lm_label=None):
        '''Processing/collecting the projected band data before the plotting function
            proj_wf = [spin,kpt,band,atom,lm] , read mcu.vasp.read.vasprun.get_projected for more details info
            
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
       
        
        # Collecting/combining the projected wfn from vasprun.xml
        if isinstance(vasprun,mcu.vasp.read.vasprun):                       # For one vasprun.xml file
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
         

        if style == 1:
            lm_shortcut = ['p','d','sp','ps','pd','dp','sd','ds','spd','sdp','psd','pds','dsp','dps']
            # Check if the lm value is appropriate
            if isinstance(lm,str):
                if lm not in lm_list and lm not in lm_shortcut:
                    raise Exception("WARNING:", lm, "is not recognizable. lm must be", lm_list, lm_shortcut)
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
                        if each_lm not in lm_list:
                            raise Exception("WARNING:", lm, "is not recognizable. lm must be one of these", lm_list)
                    else:
                        for orb in each_lm:
                            if orb  not in lm_list:
                                raise Exception("WARNING:", orb , "is not recognizable. lm must be one of these", lm_list)                        
            else:
                raise Exception("lm is not recognizable")
                
            # Compute pband
            proj_wf = proj_wf.sum(axis=2)       # Sum over the atoms --> [kpt,band,lm]
            total = proj_wf.sum(axis=2)         # Sum over the lm  --> [kpt,band]
            pband = [] 
            for each_lm in lm:
                if isinstance(each_lm,str):  
                    idx_lm = lm_list.index(each_lm)
                    proj_val = proj_wf[:,:,lm_list.index(each_lm)]/total
                else:
                    proj_val = 0
                    for orb in each_lm:
                        idx_lm = lm_list.index(orb)
                        proj_val += proj_wf[:,:,idx_lm]/total
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
                atom_list = self.vasprun.atom
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
            pband = [] 
            total = proj_wf.sum(axis=(2,3))       # Sum over the atoms --> [kpt,band,lm]
            for i in range(len(atoms)):
                idx_atom = [idx for idx in range(len(self.atom)) if atoms[i] == self.atom[idx]]
                idx_lm = [lm_list.index(lm) for lm in lms[i]] 
                proj_val_atom = 0
                proj_val = 0                
                for idx in idx_atom: proj_val_atom += proj_wf[:,:,idx,:]
                for idx in idx_lm: proj_val += proj_val_atom[:,:,idx]
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
                        if each_lm not in lm_list:
                            raise Exception("WARNING:", lm, "is not recognizable. lm must be one of these", lm_list)
                    else:
                        for orb in each_lm:
                            if orb  not in lm_list:
                                raise Exception("WARNING:", orb , "is not recognizable. lm must be one of these", lm_list) 
            else:
                raise Exception("lm is not recognizable")
                
            # Compute pband
            proj_wf = proj_wf.sum(axis=2)       # Sum over the atoms --> [kpt,band,lm]
            pband = [] 
            for each_lm in lm:                  # only two lm
                if isinstance(each_lm,str):  
                    idx_lm = lm_list.index(each_lm)
                    proj_val = proj_wf[:,:,idx_lm]
                else:
                    proj_val = 0
                    for orb in each_lm:
                        idx_lm = lm_list.index(orb)
                        proj_val += proj_wf[:,:,idx_lm]
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
           
            Attribute:
                efermi      : a Fermi level or a list of Fermi levels, it is automatically extracted frim vasprun.xml or OUTCAR
                spin        : 0  for spin unpolarized and LSORBIT = .TRUE.
                                  0 or 1 for spin polarized
                label       : label for high symmetric points, e.g. 'G-X-L-W-G',
                                  if hybridXC=True, the lavel should be a list of labels plus their coordinates
                                  
                ########################PLOTTING STYLE###################################
                style = 1   : all atoms are considered
                             lm = 's', 'py', 'pz', 'px', 'dxy', 'dyz','dz2','dxz','x2-y2' or a list of them
                                 'sp', 'pd', 'sd', 'spd'  => shortcut
                                 each color is used for each lm
                                 the marker's radius is proportional to the % of lm 
                style = 2   : considering only a list of orbitals
                             e.g. orb = ['Ni_s','C_pz']
                style = 3   : gradient map to show the character transition
                             lm = 'sp', 'pd', 'sd'
                #########################################################################
                             
                band        : the first value indicates the number of valence bands from the VBM
                              the second value indicates the number of conduction bands from the CBM
                color       : a list of strings indicating the color, following matplotlib
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

        if style == 2 and lm == 'spd' : lm = [atom+':s,p,d' for atom in self.atype]       
        if style == 3 and lm == 'spd' : lm = 'sp'   
       
        # Check if the band values are reasonable otherwise generate it
        band_idx = band
        if band_idx == None:
            idx_vbm = int(self.nelec)
            if self.lsorbit == False: idx_vbm = idx_vbm//2               # Estimation for insulator, generally wrong for metal
            first_band = int(idx_vbm - 5)
            last_band = int(idx_vbm + 4)
            if first_band < 0: first_band = 0
            if last_band > self.nbands - 1: last_band = self.nbands - 1
            band_idx = [first_band, last_band]
        else:
            assert band_idx[0] <= band_idx[1]                    # from band[0] to band[1]
            if band_idx[0] < 1: band_idx[0] = 1     # index used in OUTCAR, will be shifted to start at zero
            if band_idx[1] > self.nbands: band_idx[1] = self.nbands              # Cannot larger than the number of bands
            band_idx[0] = band_idx[0] -1
            band_idx[1] = band_idx[1] -1     
        
        band, path, sym_kpoint_coor, label, conventional = self._generate_band(self.vasprun, efermi, spin, label)  
        pband = self._generate_pband(self.vasprun, spin, style, lm)
        
        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        yrange = (-50,50)
        
        # Customization:
        border = 1.08
        
        # Plot the high symmetric kpoint grid
        if conventional == True or label != None:
            for kpt in range(sym_kpoint_coor.shape[0]):
                ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color=band_color[1],linewidth=1.0)

        if label != None and xlim == None:
            nkpts = len(label)
            assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
            for kpt in range(nkpts):   
                point = label[kpt]
                if point == 'G': point = r'$\Gamma$'
                ax.text(sym_kpoint_coor[kpt]/path.max()+0.015, -0.065, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                        color='black', fontsize=fontsize)     
            
        # Plot bands            
        ax.plot([0, path.max()], [0,0], color=band_color[2], linewidth=1.0, dashes=[6,3])
        for ith in range(band.shape[1]):
            ax.plot(path.T, band[:,ith], color=band_color[0],linewidth=1.0)    
             
        # Plot pbands 
        color_list = ['r','g','b','y','m','c']
        path = np.array(path).flatten() 
        if style == 1 or style == 2:
            pband = 200 * scale * np.power(pband,2)     # The radius of the marker ~ the percent 
            
            # Color
            if color == None: 
                color = color_list
            else:
                assert isinstance(color,list) or isinstance(color,str)
                if isinstance(color,str): color = [color]
                
            if facecolor == None: 
                fcolors = color
            elif isinstance(facecolor,list):
                assert len(facecolor) == len(pband)
                fcolors = []
                for i in range(len(facecolor)):
                    assert facecolor[i] == True or facecolor[i] == False
                    if facecolor[i] == True: fcolors.append(color[i]) 
                    if facecolor[i] == False: fcolors.append('none') 
            elif facecolor == 'none':
                fcolors = ['none']*len(pband)

            # Marker
            if marker == None: 
                marker = ['o']*len(pband)
            else:
                assert isinstance(marker,list) or isinstance(legend,str)
                if isinstance(marker,str): marker = [marker]
                assert len(marker) == len(pband)                
            
            # legend    
            if legend != None:
                legends = []   
                assert isinstance(legend,list) or isinstance(legend,str)
                if isinstance(legend,str): legend = [legend]
                assert len(legend) == len(pband)
                
            # Actual plotting
            for lm in range(len(pband)):
                for ith in range(band_idx[0],band_idx[1]):
                    ax.scatter(path, band[:,ith], s=pband[lm][:,ith], facecolors=fcolors[lm], edgecolors=color[lm], alpha=alpha, marker=marker[lm])
                ith = band_idx[1]
                if legend == None:
                    ax.scatter(path, band[:,ith], s=pband[lm][:,ith], facecolors=fcolors[lm], edgecolors=color[lm], alpha=alpha, marker=marker[lm])
                else:
                    ax.scatter(path, band[:,ith], s=pband[lm][:,ith], facecolors=fcolors[lm], edgecolors=color[lm], alpha=alpha, marker=marker[lm],label=legend[lm])                    
                
            if legend != None: 
                lgnd = ax.legend(loc=loc, numpoints=1, fontsize=fontsize)
                for i in range(len(pband)): lgnd.legendHandles[i]._sizes = [legend_size*60]
                
        elif style == 3:
            path = np.array(path).flatten()
            if marker == None: 
                marker = 'o'
            else:
                assert isinstance(marker,str)
            for ith in range(band_idx[0],band_idx[1]+1):
                plt.scatter(path, band[:,ith], c=pband[:,ith], s=50*scale, vmin=0.0, vmax=1., cmap=cmap, marker=marker, edgecolor=edgecolor) 
            cbar = plt.colorbar(ticks=[])
            cbar.outline.set_linewidth(border)
        
        # Graph adjustments             
        ax.tick_params(labelsize=fontsize, width=border)
        ax.spines['top'].set_linewidth(border)
        ax.spines['right'].set_linewidth(border)
        ax.spines['bottom'].set_linewidth(border)
        ax.spines['left'].set_linewidth(border)
        if xlim == None:
            plt.xlim([0,path.max()])
            plt.xticks([])
            plt.xlabel('k', size=fontsize+4)
        else:
            plt.xlim(xlim)
            plt.xlabel('k ' + r'($\AA^{-1}$)', size=fontsize+4)
        ax.xaxis.set_label_coords(0.5, -0.08) 
        plt.ylabel('Energy (eV)', size=fontsize+4)        
        plt.ylim(ylim)
        plt.tight_layout()
        if save == True: 
            fig.savefig(figname+'.'+format, dpi=dpi, format=format)      
        else:
            plt.show() 
            

    def _generate_dos(self, vasprun, efermi=None, spin=0, lm=None):
        '''Processing/collecting the DOS data before the plotting function
            Note: unlike plot_band function, only one vasprun.xml is used. Combining vasprun.xml for DOS sounds a bit weird
            and unececessary
            
            spin            : spin of DOS.
                              For LSORBIT == True: spin = 0,1,2,3
                              For ISPIN = 2      : spin = 0,1
                              
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        
        # Get the fermi level
        pdos_exist = False
        lm_list = vasprun.lm
        
        vasprun.get_dos()
        tdos = vasprun.tdos[spin,:,:2]
        if vasprun.pdos_exist == True: 
            pdos_data = vasprun.pdos[:,spin,:,:]         # [atom,energy,lm]
            pdos_exist = True
        else:
            pdos = None
            
             
        # Compute pDOS
        if pdos_exist == True:
        
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
                atom_list = self.vasprun.atom
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
            
            # Recompute tDOS from pdos_data, the one provided in the vasprun does not neccessarily equal to this tDOS
            # however, this one is consistent with the pDOS below
            temp = pdos_data[:,:,1:].sum(axis=0)
            tdos = np.empty([temp.shape[0],2])
            tdos[:,0] = pdos_data[0,:,0]          
            tdos[:,1] = temp.sum(axis=1)  
            
            # Compute pDOS
            pdos = [] 
            for i in range(len(atoms)):
                idx_atom = [idx for idx in range(len(self.atom)) if atoms[i] == self.atom[idx]]
                idx_lm = [lm_list.index(lm) for lm in lms[i]] 
                proj_val_atom = 0
                proj_val = 0
                for idx in idx_atom: proj_val_atom += pdos_data[idx,:,1:]        # Sum over all atoms
                for idx in idx_lm: proj_val += proj_val_atom[:,idx]
                pdos.append(proj_val)
            pdos = np.asarray(pdos).T
            
        # Shift the energy 
        tdos[:,0] = tdos[:,0] - efermi
                
        return tdos, pdos_exist, pdos
        
    def plot_dos(self, vasprun=None, style=1, efermi=None, spin=0, lm=None, color=None,
                    legend=None, loc="upper right", fill=True, alpha=0.2,
                    save=False, figname='DOS', figsize=(6,3), elim=(-6,6), yscale=1.1, fontsize=18, dpi=600, format='png'):
        '''Plot projected band structure
            For multiple vasprun.xml, user can choose one of them to plot the DOS. Default: the first vasprun.xml

            Attribute:
            style           : 1 (standard plot) or 2 (vertital plot)
            
            spin            : spin of DOS.
                              For LSORBIT == True: spin = 0,1,2,3
                              For ISPIN = 2      : spin = 0,1
                              
            lm              : string or a list of string, e.g. 'Ni:s' or ['Ni:s','C:s,px,pz']
        '''
        
        if vasprun == None: 
            if isinstance(self.vasprun,mcu.vasp.read.vasprun): 
                vasprun = self.vasprun
                if efermi == None: efermi = self.efermi
            if isinstance(self.vasprun,list): 
                vasprun = self.vasprun[0]  
                if efermi == None: efermi = self.efermi[0]
        else:
            assert isinstance(vasprun,mcu.vasp.read.vasprun)
            
        if lm == None: 
            lm = [atom+':s,p,d' for atom in self.atype]  
            legend = lm
        elif lm != None and legend == None:
            legend = lm
              
        if spin == 'updown':
            if self.ispin != 2: raise Exception('ISPIN must be 2 for the up-down DOS plotting')
            tdos0, pdos_exist, pdos0 = self._generate_dos(vasprun, efermi=efermi, spin=0, lm=lm)
            tdos1, pdos_exist, pdos1 = self._generate_dos(vasprun, efermi=efermi, spin=1, lm=lm)
            tdos1 = -tdos1
            pdos1 = -pdos1
            if figsize == (6,3) and style==1 : figsize = (6,5)
            if figsize == (6,3) and style==2 : figsize = (5,6)
        else:
            if figsize == (6,3) and style==2 : figsize = (3,6)
            tdos0, pdos_exist, pdos0 = self._generate_dos(vasprun, efermi=efermi, spin=spin, lm=lm)
        
        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        color_list = ['k','r','g','b','y','m','c']
        if color == None: color = color_list
        
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        yrange = (-50,50)
        
        # Plot DOS    
        if style == 1:
            ax.plot(tdos0[:,0], tdos0[:,1], color=color[0],linewidth=1.1,label='TDOS')
            if pdos_exist == True:
                for orb in range(pdos0.shape[1]): 
                    ax.plot(tdos0[:,0], pdos0[:,orb], color=color[orb+1],linewidth=1.0,label=legend[orb])
                    if fill == True: ax.fill(tdos0[:,0], pdos0[:,orb], color=color[orb+1], alpha=alpha)
                
            if spin == 'updown':
                ax.plot(tdos0[:,0], tdos1[:,1], color=color[0],linewidth=1.1)
                if pdos_exist == True:
                    for orb in range(pdos1.shape[1]): 
                        ax.plot(tdos0[:,0], pdos1[:,orb], color=color[orb+1],linewidth=1.0)
                        if fill == True: ax.fill(tdos0[:,0], pdos1[:,orb], color=color[orb+1], alpha=alpha)
        
            # Graph adjustments 
            plt.xlabel('Energy (eV)', size=fontsize+4)   
            plt.ylabel('DOS', size=fontsize+4)
            if spin == 'updown':
                plt.ylim([tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale])  
                ax.plot([0,0], [tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
                ax.plot([tdos0[:,0].min()*yscale,tdos0[:,0].max()*yscale], [0,0], color=color[0], linewidth=1.0, alpha=alpha) 
            else:
                plt.ylim([0, tdos0[:,1].max()*yscale])
                ax.plot([0,0], [0, tdos0[:,1].max()*yscale], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
            plt.xlim(elim)
            plt.yticks([])
            
        elif style == 2:
            ax.plot(tdos0[:,1], tdos0[:,0], color=color[0],linewidth=1.1,label='TDOS')
            if pdos_exist == True:
                for orb in range(pdos0.shape[1]): 
                    ax.plot(pdos0[:,orb], tdos0[:,0], color=color[orb+1],linewidth=1.0,label=legend[orb])
                    if fill == True: ax.fill(pdos0[:,orb], tdos0[:,0], color=color[orb+1], alpha=alpha)
                
            if spin == 'updown':
                ax.plot(tdos1[:,1], tdos0[:,0], color=color[0],linewidth=1.1)
                if pdos_exist == True:
                    for orb in range(pdos1.shape[1]): 
                        ax.plot(pdos1[:,orb], tdos0[:,0], color=color[orb+1],linewidth=1.0)
                        if fill == True: ax.fill(pdos1[:,orb], tdos0[:,0], color=color[orb+1], alpha=alpha)

        
            # Graph adjustments 
            plt.xlabel('DOS', size=fontsize+4)   
            plt.ylabel('Energy (eV)', size=fontsize+4)
            if spin == 'updown':
                plt.xlim([tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale])  
                ax.plot([tdos1[:,1].min()*yscale, tdos0[:,1].max()*yscale], [0,0], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
                ax.plot([0,0], [tdos0[:,0].min()*yscale,tdos0[:,0].max()*yscale], color=color[0], linewidth=1.0, alpha=alpha) 
            else:
                plt.xlim([0,tdos0[:,1].max()*yscale])
                ax.plot([0, tdos0[:,1].max()*yscale], [0,0], color=color[0], linewidth=1.0, dashes=[6,3], alpha=alpha) 
            plt.ylim(elim)
            plt.xticks([])
            
            
        # Legend
        lgnd = ax.legend(loc=loc, numpoints=1, fontsize=fontsize)
                
        # Graph adjustments 
        border = 1.08        
        ax.tick_params(labelsize=fontsize, width=border)
        ax.spines['top'].set_linewidth(border)
        ax.spines['right'].set_linewidth(border)
        ax.spines['bottom'].set_linewidth(border)
        ax.spines['left'].set_linewidth(border)
        plt.tight_layout()
        if save == True: 
            fig.savefig(figname+'.'+format, dpi=dpi, format=format)      
        else:
            plt.show() 
            
    def _generate_spin(self, vasprun, lm=None):
        '''Processing/collecting the spin texture data before the plotting function
            proj_wf = [spin,kpt,band,atom,lm] , read mcu.vasp.read.vasprun.get_projected for more details info
            
            lm          : ['Ni:s','C:pz']
        '''      
        
        # Collecting/combining the projected wfn from vasprun.xml
        if isinstance(vasprun,mcu.vasp.read.vasprun):                       # For one vasprun.xml file
            vasprun.get_projected()
            lm_list = vasprun.lm 
            proj_wfs = []
            for spin in range(1,4):
                proj_wf = vasprun.proj_wf[spin] 
                #proj_wf = utils.rm_redundant_band(self.kpts, proj_wf)[1]          # remove redundant 
                weight = vasprun.kpoints['weights']
                nonzero = np.count_nonzero(weight)
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
            atom_list = self.vasprun.atom
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
            
        # Compute spin texture
        spin_text = [] 
        for spin in range(3):
            spin_ = []
            for i in range(len(atoms)):
                idx_atom = [idx for idx in range(len(self.atom)) if atoms[i] == self.atom[idx]]
                idx_lm = [lm_list.index(lm) for lm in lms[i]] 
                proj_val_atom = 0
                proj_val = 0                
                for idx in idx_atom: proj_val_atom += proj_wfs[spin][:,:,idx,:]
                for idx in idx_lm: proj_val += proj_val_atom[:,:,idx]
                spin_.append(proj_val)
            spin_ = np.asarray(spin_).sum(axis=0)
            spin_text.append(spin_)
        
        return np.asarray(spin_text)

    def plot_spin(self, style=1, lm=None, band=None, cmap='bwr',
                    save=False, figname='spin_texture', figsize=(6,6), xlim=None, ylim=None, fontsize=18, dpi=600, format='png'):
        '''Plot spin texture
           
            Attribute:
                efermi      : a Fermi level or a list of Fermi levels, it is automatically extracted frim vasprun.xml or OUTCAR
                lm          : 'Ni:s' or ['Ni:s','C:pz']
                
                band        : index of the band
                
                color       : a list of strings indicating the color, following matplotlib
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

        if lm == None : lm = [atom+':s,p,d' for atom in self.atype]         
        
        # Check if the band values are reasonable otherwise generate it
        if band == None:
            idx_vbm = int(self.nelec)
            if self.lsorbit == False: idx_vbm = idx_vbm//2
            band = idx_vbm - 1
        else:
            assert (band <= self.nbands) and (band >= 1)  
            band = band - 1

        spin_text = self._generate_spin(self.vasprun, lm=lm)[:,:,band]

        # Get X, Y
        kpoint = read.KPOINTS()
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
   
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        
        if style == 2:
            my_cm = plt.cm.bwr
            # colors = ['#0080ff', '#ffffff', '#ff0066']
            # my_cm = LinearSegmentedColormap.from_list('my_list', colors, N=100)
            plt.contour(X, Y, C, linewidths=0.01)
            plt.contourf(X, Y, C, vmin=-1.0, vmax=1.0, cmap='bwr')
            cbar = plt.colorbar(ticks=[-1.0,0.0,1.0])
            cbar.ax.tick_params(labelsize=fontsize) 
            ax.quiver(X, Y, U, V, color='k')
        elif style == 1:
            ax.quiver(X, Y, U, V, C, angles='xy')
        
        # Graph adjustments             
        ax.tick_params(labelsize=fontsize)
        if xlim == None: xlim = [X.min(),X.max()]
        if ylim == None: ylim = [Y.min(),Y.max()]
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
        ax.set_xticks([xlim[0],0,xlim[1]])
        ax.set_yticks([ylim[0],0,ylim[1]])
        plt.tight_layout()
        if save == True: 
            fig.savefig('Band.'+format,dpi=dpi,format=format)      
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
        
        assert isinstance(band_color,list)
        assert len(band_color) == 3

        band_idx = band
        if band_idx == None:
            idx_vbm = int(self.nelec)
            if self.lsorbit == False: idx_vbm = idx_vbm//2               # Estimation for insulator, generally wrong for metal
            band_idx = [idx_vbm-1, idx_vbm]
        else:
            assert band_idx[0] <= band_idx[1]                    # from band[0] to band[1]
            if band_idx[0] < 1: band_idx[0] = 1     # index used in OUTCAR, will be shifted to start at zero
            if band_idx[1] > self.nbands: band_idx[1] = self.nbands              # Cannot larger than the number of bands
            band_idx[0] = band_idx[0] -1
            band_idx[1] = band_idx[1] -1   
            
        band, path, sym_kpoint_coor, label, conventional = self._generate_band(self.vasprun, efermi, spin, label=None)  
        
        # Get X, Y
        kpoint = read.KPOINTS()
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
            
class LOCPOT:
    def __init__(self, locpot='LOCPOT'):
        '''Get LOCPOT file and return a LOCPOT object '''
        self.locpot = read.LOCPOT(locpot)
        
    def get_2D_average(self, axis='z'):
        return self.locpot.get_2D_average(axis)
        
    def get_vacumm(self, pot=None, axis='z', error=0.01):
        '''Get the electrostatic potential at vacuum
        '''
        
        if not isinstance(pot,np.ndarray): 
            pot = self.get_2D_average(axis)
        lower_bound = pot[1].max()- 2*error
        idx = pot[1] > lower_bound
        pot_in_window = pot[1,idx]
        e_vacuum = np.average(pot[1,idx])
    
        return e_vacuum
        
    def plot(self, axis='z', error=0.01, color=['r', '#737373'], ylim=None, save=False, figname='elecpot', figsize=(8,6), fontsize=18, dpi=600, format='png'):
        '''Function to plot the inplane average potential to check the convegence
        '''
        
        pot = self.get_2D_average(axis)
        e_vacuum = self.get_vacumm(pot, error=error)
        
        ##----------------------------------------------------------
        ##Plotting:        
        ##----------------------------------------------------------
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        ax.plot(pot[0], pot[1], color=color[0],linewidth=1.1,label='TDOS')
        ax.plot([pot[0].min(), pot[0].max()], [e_vacuum,e_vacuum], color=color[1], linewidth=1.1, dashes=[6,3])
        if ylim == None: ylim = (round(pot[1].min()) - 10.0,round(pot[1].max()) + 10.0)
        y_evacuum = (e_vacuum-ylim[0])/(ylim[1]-ylim[0]) 
        label = 'Vacuum = ' + str(round(e_vacuum,2)) + r' $\pm$ ' + str(round(error,2)) + ' eV' 
        ax.text(0.02, y_evacuum, label  , verticalalignment='bottom', horizontalalignment='left',transform=ax.transAxes, color='black', fontsize=fontsize+4)  
 
        # Graph adjustments
        border = 1.08
        ax.tick_params(labelsize=fontsize, width=border)
        ax.spines['top'].set_linewidth(border)
        ax.spines['right'].set_linewidth(border)
        ax.spines['bottom'].set_linewidth(border)
        ax.spines['left'].set_linewidth(border)
        plt.xlabel(axis + r' ($\AA$)', size=fontsize+4)
        ax.xaxis.set_label_coords(0.5, -0.09) 
        plt.ylabel('Electrostatic potential (V)', size=fontsize+4)        
        plt.ylim(ylim) 
        plt.xlim([0,pot[0].max()])
        plt.tight_layout()
        if save == True: 
            fig.savefig(figname+'.'+format, dpi=dpi, format=format)      
        else:
            plt.show() 
            
class POSCAR:
    def __init__(self, poscar='POSCAR'):
        '''Get POSCAR file and return a POSCAR object '''
        self.poscar = read.POSCAR(poscar)
        self.cell_recip = self.poscar.cell[1]
        
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
                

        
    