#!/usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import mcu
from mcu.vasp import utils, vasprun

class VASP:
    def __init__(self, path='./', vaspruns='vasprun'):
        '''
            path        : the project directory
            vaspruns    : a str or a list of string as names for *.xml files
        '''
        
        # Create vasprun object(s)
        if isinstance(vaspruns, str):
            self.vasprun = vasprun.read(path + '/' + vaspruns + '.xml')
        elif isinstance(vaspruns, list):
            self.vasprun = []
            for xml in vaspruns:
                xml_file = path + '/' + xml + '.xml'
                if not utils.check_exist(xml_file):
                    print('Cannot find:', xml_file)
                    break
                self.vasprun.append(vasprun.read(xml_file))
            
        else:
            print('Provide a string or a list of names for *.xml file')
    
        
    def get_bandgap(self, efermi=0):
        '''Get the bandgap'''

        if isinstance(self.vasprun,mcu.vasp.vasprun.read):
            self.vasprun.get_band()
            self.band = self.vasprun.band[:,:,:,0]
            self.co_occ = self.vasprun.band[:,:,:,1] 
            self.co_occ_ = self.co_occ > 0.5       
            electronic = self.vasprun.parameters['electronic']
            ispin = electronic.spin['ISPIN']
            self.kpts = self.vasprun.kpoints['kpointlist']
            nkpts = self.kpts.shape[0]
            
        elif isinstance(self.vasprun,list):
            if efermi == 0:
                efermi = [0.0]*len(self.vasprun)
            elif not isinstance(efermi,list): 
                raise Exception('Please provide a list of Fermi level for each *.xml')
                
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
                
                
    def plot_band(self, efermi=0, spin=0, label=None, hybridXC=False, save=False,
                    figsize=(5,8), ylim=[-6,6], fontsize=18, dpi=600, format='png'):
        '''Plot band structure
           
            Attribute:
                efermi          : a Fermi level or a list of Fermi levels
                spin            : either 0 (close shell) or 1 (spin-polarized)
                label           : label for high symmetric points, e.g. G-X-L-W-G 
        
        
        '''
        
        # matplotlib is required
        import matplotlib
        import matplotlib.pyplot as plt
        
        coor_highkpt = None
        if label == None:
            raise Exception('Please provide label for high symmetric k-points')
        elif hybridXC == False:
            label = label.split('-')
        elif hybridXC == True:
            if not isinstance(label,list): raise Exception('Please provide label and coordinates for high symmetric k-points')
            temp = []
            coor_kpts = [] 
            for kpt in label:
                temp.append(kpt[0])
                coor_kpts.append(kpt[1:])
            label = temp       
            coor_kpts = np.asarray(coor_kpts)
                
        nkpts = len(label)             
        band = None
        path = None
        
        if hybridXC == False and isinstance(self.vasprun,mcu.vasp.vasprun.read): # For conventional band structure calculation 
            self.vasprun.get_band()
            band = self.vasprun.band[spin][:,:,0]
            kpts = self.vasprun.kpoints['kpointlist']
            kpts, band = utils.rm_redundant_band(kpts, band)
            
            # Find absolute kpts and shift the band
            b = self.vasprun.cell_final[1]          # Get the reciprocal lattice
            abs_kpts = kpts.dot(b)                  # From fractional to absolute
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            path = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())
            band = band - efermi               # Efermi is set at 0 eV
            
            highsym_kpt = self.vasprun.kpoints['points']
            if nkpts != highsym_kpt.shape[0]: 
                raise Exception('The no. of kpoint label provided does not match the no. of path in the calculation')
            else:
                coor_highkpt = [0.0]
                for kpt in range(nkpts-2):
                    idx = ((path.shape[1] + nkpts - 2)//(nkpts-1) - 1) * (kpt+1)
                    coor = path[0,idx]         
                    coor_highkpt.append(coor)
                coor_highkpt.append(1.0*path.max())   
                coor_highkpt = np.asarray(coor_highkpt)
                     
        elif hybridXC == True:
            if isinstance(self.vasprun,mcu.vasp.vasprun.read):
                vasprun = self.vasprun
                vasprun.get_band()
                band = vasprun.band[spin][:,:,0]
                kpts = vasprun.kpoints['kpointlist']
                weight = vasprun.kpoints['weights']
                nonzero = np.count_nonzero(weight)
                kpts, band = kpts[nonzero:], band[nonzero:]
                band = band - efermi
            elif isinstance(self.vasprun,list):
                if efermi == 0:
                    efermi = [0.0]*len(self.vasprun)
                elif not isinstance(efermi,list): 
                    raise Exception('Please provide a list of Fermi level for each *.xml')
                    
                electronic = self.vasprun[0].parameters['electronic']
                nbands = electronic.general['NBANDS']
                bands = np.zeros([1,nbands])
                kptss = np.zeros([1,3])
                for i, vasprun in enumerate(self.vasprun):
                    vasprun.get_band()
                    band = vasprun.band[spin][:,:,0]
                    kpts = vasprun.kpoints['kpointlist']
                    weight = vasprun.kpoints['weights']
                    nonzero = np.count_nonzero(weight)
                    kpts, band = kpts[nonzero:], band[nonzero:]
                    band = band - efermi[i]
                    bands = np.vstack([bands,band])
                    kptss = np.vstack([kptss,kpts])
                    
                kpts, band = kptss[1:], bands[1:]
                vasprun = self.vasprun[0]
                                
            # Find absolute kpts
            b = vasprun.cell_final[1]          # Get the reciprocal lattice
            abs_kpts = kpts.dot(b)                  # From fractional to absolute
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            path = np.matrix(np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum())

            # Find absolute coordinates for high symmetric kpoints     
            abs_kpts = coor_kpts.dot(b)   
            temp_kpts = np.empty_like(abs_kpts)
            temp_kpts[0] = abs_kpts[0]
            temp_kpts[1:] = abs_kpts[:-1] 
            coor_highkpt = np.sqrt(((temp_kpts - abs_kpts)**2).sum(axis=1)).cumsum() 

        # Plotting:            
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
        yrange = (-50,50)
        for kpt in range(nkpts-2):
            ax.plot([coor_highkpt[kpt+1]]*2,yrange,color="#808080",linewidth=1.0)

            
        for kpt in range(nkpts):   
            point = label[kpt]
            if point == 'G': point = r'$\Gamma$'
            ax.text(coor_highkpt[kpt]/path.max()+0.02, -0.05, point, verticalalignment='bottom', horizontalalignment='right',transform=ax.transAxes,
                    color='black', fontsize=fontsize)     
            
        ax.plot([0,path.max()],[0,0],color="#808080",linewidth=1.0, dashes=[6,3])
        for ith in range(band.shape[1]):
            ax.plot(path.T,band[:,ith],color="#007acc",linewidth=1.0)    
                
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
            
        
    def plot_proband(self):
        pass        
        
    def plot_dos(self):
        pass       
      
        
        

