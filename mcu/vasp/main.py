#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import mcu
from mcu.vasp import utils, read

class VASP:
    def __init__(self, path='./', vaspruns='vasprun', outcars=None):
        '''
            path        : the project directory
            vaspruns    : a str or a list of string as names for *.xml files
            outcars     : a str or a list of string as names for OUTCAR files
        '''
        
        # Create vasprun object(s)
        if isinstance(vaspruns, str):                   # For one vasprun.xml file    
            self.vasprun = read.vasprun(path + '/' + vaspruns + '.xml')
            self.useOUTCAR = False
            if outcars != None: 
                self.outcar = read.OUTCAR(outcars)
                self.useOUTCAR = True
        elif isinstance(vaspruns, list):                # For multiple vasprun.xml file
            self.vasprun = []
            for xml in vaspruns:
                xml_file = path + '/' + xml + '.xml'
                if not utils.check_exist(xml_file):
                    print('Cannot find:', xml_file)
                    break
                self.vasprun.append(read.vasprun(xml_file))
                
            if outcars != None: 
                assert isinstance(outcars, list)               
                assert len(outcars) == len(vaspruns)
                self.outcar = []
                for outcar in outcars:
                    outcar_file = path + '/' + outcar
                    if not utils.check_exist(outcar_file):
                        print('Cannot find:', outcar_file)
                        break
                    self.outcar.append(read.OUTCAR(outcar_file))  
                    self.useOUTCAR = True                    
        else:
            print('Provide a string or a list of names for *.xml file')
    
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
        if efermi == None: 
            self.get_efermi()
            efermi = self.efermi
            
        if isinstance(self.vasprun,mcu.vasp.read.vasprun):              # For one vasprun.xml file
            assert isinstance(efermi,float) 
            self.vasprun.get_band()
            self.band = self.vasprun.band[:,:,:,0]
            self.co_occ = self.vasprun.band[:,:,:,1] 
            self.co_occ_ = self.co_occ > 0.5       
            electronic = self.vasprun.parameters['electronic']
            ispin = electronic.spin['ISPIN']
            self.kpts = self.vasprun.kpoints['kpointlist']
            nkpts = self.kpts.shape[0]
        elif isinstance(self.vasprun,list):                             # For multiple vasprun.xml file
            assert isinstance(efermi,list)
            for i in range(len(self.vasprun)): 
                assert isinstance(efermi[i],float)
                
            electronic = self.vasprun[0].parameters['electronic']
            nbands = electronic.general['NBANDS']
            ispin = electronic.spin['ISPIN']
            
            band_spin = [] 
            co_occ_spin1 = []
            co_occ_spin2 = []            
            for spin in range(ispin):
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
            
        bandedge = np.zeros([ispin,nkpts,2,2])
        self.bandgap = []
        for spin in range(ispin):
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
        if efermi == None: 
            self.get_efermi()
            efermi = self.efermi        
        
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
            b = vasprun.cell_final[1]          # Get the reciprocal lattice
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
                vasprun = vasprun
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
            b = vasprun.cell_final[1]          # Get the reciprocal lattice
            abs_kpts = kpts.dot(b)                  # From fractional to absolute
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
                
    def plot_band(self, efermi=None, spin=0, label=None, save=False,
                    figsize=(5,8), ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : either 0 (spin unpolarized) or 1 (spin polarized)
                label           : label for high symmetric points, e.g. G-X-L-W-G,
                                  if hybridXC=True, the lavel should be a list of labels plus their coordinates
        '''
        
        # matplotlib is required
        import matplotlib
        import matplotlib.pyplot as plt
            
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
                ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color="#808080",linewidth=1.0)

        if label != None:
            nkpts = len(label)
            assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
            for kpt in range(nkpts):   
                point = label[kpt]
                if point == 'G': point = r'$\Gamma$'
                ax.text(sym_kpoint_coor[kpt]/path.max()+0.02, -0.05, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                        color='black', fontsize=fontsize)     
            
        # Plot bands            
        ax.plot([0,path.max()],[0,0],color="#808080",linewidth=1.0, dashes=[6,3])
        for ith in range(band.shape[1]):
            ax.plot(path.T,band[:,ith],color="#007acc",linewidth=1.0)    
             
        # Graph adjustments             
        ax.tick_params(labelsize=fontsize)
        plt.xlabel('k', size=fontsize+4)
        ax.xaxis.set_label_coords(0.5, -0.08) 
        plt.ylabel('Energy (eV)', size=fontsize+4)        
        plt.ylim(ylim)
        plt.xlim([0,path.max()])
        plt.xticks([])
        plt.tight_layout()
        if save == True: 
            fig.savefig('Band.'+format,dpi=dpi,format=format)      
        else:
            plt.show()            
            
        
    def _generate_pband(self, vasprun):
        '''Processing/collecting the projected band data before the plotting function
        '''        
        if isinstance(vasprun,mcu.vasp.read.vasprun):                       # For one vasprun.xml file
            vasprun.get_projected()
            pband = vasprun.pro_wf[0]       [kpt,band,atom,l]
        elif isinstance(vasprun,list):                                      # For multiple vasprun.xml file
            pass
            # for i, run in enumerate(vasprun):
                # run.get_projected()
                # pband = run.pro_wf[0][:,:,0]
                # kpts = run.kpoints['kpointlist']
                # weight = run.kpoints['weights']
                # nonzero = np.count_nonzero(weight)
                # kpts, band = kpts[nonzero:], band[nonzero:]
                    
                    
                    
    def plot_pband(self, efermi=None, spin=0, label=None, save=False,
                    figsize=(5,8), ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot projected band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : either 0 (spin unpolarized) or 1 (spin polarized)
                label           : label for high symmetric points, e.g. G-X-L-W-G,
                                  if hybridXC=True, the lavel should be a list of labels plus their coordinates
            
        '''
                    
        # matplotlib is required
        import matplotlib
        import matplotlib.pyplot as plt
            
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
                ax.plot([sym_kpoint_coor[kpt]]*2,yrange,color="#808080",linewidth=1.0)

        if label != None:
            nkpts = len(label)
            assert nkpts == sym_kpoint_coor.shape[0]        # The numbers of label should be match with the # of coordiantes provided
            for kpt in range(nkpts):   
                point = label[kpt]
                if point == 'G': point = r'$\Gamma$'
                ax.text(sym_kpoint_coor[kpt]/path.max()+0.02, -0.05, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                        color='black', fontsize=fontsize)     
            
        # Plot bands            
        ax.plot([0,path.max()],[0,0],color="#808080",linewidth=1.0, dashes=[6,3])
        for ith in range(band.shape[1]):
            ax.plot(path.T,band[:,ith],color="#007acc",linewidth=1.0)    
             
        # Graph adjustments             
        ax.tick_params(labelsize=fontsize)
        plt.xlabel('k', size=fontsize+4)
        ax.xaxis.set_label_coords(0.5, -0.08) 
        plt.ylabel('Energy (eV)', size=fontsize+4)        
        plt.ylim(ylim)
        plt.xlim([0,path.max()])
        plt.xticks([])
        plt.tight_layout()
        if save == True: 
            fig.savefig('Band.'+format,dpi=dpi,format=format)      
        else:
            plt.show() 
            
        
    def plot_dos(self):
        pass       
      
        
        

