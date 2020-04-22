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
from ..utils.misc import check_exist
           
        
def read_BAND(filename):
    '''Read seedname_dat.BAND/BAND.DAT file'''
    
    if not check_exist(filename):
        if check_exist("BAND.DAT"):
            print('Cannot find : ' + filename, '. BAND.DAT is used instead')
            filename = "BAND.DAT"
        else:
            assert False, 'Cannot find the BAND file. Check the path: ' + filename
    
    with open(filename, "r") as data_file:
        data = data_file.readlines()
        nkpt, nband, nspin = np.int64([data[0].split()[i] for i in [2,4,6]])
        npath = int(data[1].split()[2])
        start = npath + 3
        block_length = 19 + 2*npath + nkpt
        band = []
        for spin in range(nspin):
            block = data[(start + block_length*spin):(start + block_length*(spin + 1))]
            symm_k_coor = np.float64([block[14 + 2*i].split()[-1] for i in range(npath + 1)])
            band_block = block[(18 + 2*npath):-1]
            band_spin = np.float64("\n".join(band_block).split())
            band_spin = band_spin.reshape(nkpt, nband+1)
            band.append(band_spin[:,1:])
            kpath_abs = band_spin[:,0] 
            efermi = np.float64(block[-1].split()[-1])

        band = np.float64(band)
            
    return symm_k_coor, kpath_abs, band, efermi

def read_DOSS(filename):
    '''Read seedname_dat.DOSS/DOSS.DAT file'''
    
    if not check_exist(filename):
        if check_exist("DOSS.DAT"):
            print('Cannot find : ' + filename, '. DOSS.DAT is used instead')
            filename = "DOSS.DAT"
        else:
            assert False, 'Cannot find the DOSS file. Check the path: ' + filename

    with open(filename, "r") as data_file:
        data = data_file.readlines()
        nepts, nproj, nspin = np.int64([data[0].split()[i] for i in [2,4,6]]) 
        block_length = nepts + 2
        
        dos = []
        for spin in range(nspin):
            block = data[(3 + block_length*spin):(3 + block_length*(spin + 1))]
            dos_block = block[1:-1]
            dos_spin = np.float64("".join(dos_block).split())
            dos_spin = dos_spin.reshape(nepts, nproj + 1)
            dos.append(dos_spin)
            epts = dos_spin[:,0] 
            efermi = np.float64(block[-1].split()[-1])
            
        dos = np.float64(dos)
    
    return epts, dos, efermi
        
def read_f25(filename=None):
    '''Read seedname_dat.f25/fort.25 file'''
    if filename is None: filename = "fort.25"
    
    if not check_exist(filename):
        print('Cannot find : ' + filename, '. Searching for fort.25 ...')
        assert check_exist("fort.25"), "Cannot find the fort.25 file"
        print('Found fort.25. Band structure is extracted from fort.25')
        filename = "fort.25"
        
    with open(filename, "r") as data_file:
        data = data_file.read().split('\n')
        
        jump = 0
        out = []
        EOS = False
        while EOS == False:
            first_line = data[jump]
            ihferm = int(first_line[3])
            type = first_line[4:8]
            nband, nkp = int(first_line[8:13]), int(first_line[13:18])
            dk, efermi = float(first_line[-24:-12]), float(first_line[-12:])
            length = int(np.ceil(nband * nkp / 6))
            data_block = data[(jump+3):(jump+3+length)]
            data_block = "".join(data_block)
            eigenvals = np.float64(list(map(''.join, zip(*[iter(data_block)]*12))))
            out.append([ihferm, type, nband, nkp, dk, efermi, eigenvals])
            jump += length + 3
            
            if data[jump] is '': EOS = True
    return out
        