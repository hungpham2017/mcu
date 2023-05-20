#!/usr/bin/env python -u 
'''
pDMET: Density Matrix Embedding theory for Periodic Systems
Copyright (C) 2018 Hung Q. Pham. All Rights Reserved.
A few functions in pDMET are modifed from QC-DMET Copyright (C) 2015 Sebastian Wouters

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

import os, sys, datetime
import numpy as np        


def check_exist(file):
    '''Check if a file exists in the running directory '''        
    exist = os.path.exists(file)
    return exist
    
def print_msg(msg=None, *kargs):
    if msg == None:
        print() 
    else:
        print(msg, *kargs)  
    sys.stdout.flush()
    
def date():  
    return datetime.datetime.now().strftime("%Y/%m/%d - %H:%M:%S")    
    
def unique(array):
    array = np.asarray(array)
    unique_array = []
    counts = []
    for element in range(array.shape[0]):
        if array[element] not in unique_array: 
            unique_array.append(array[element])
            count = 1
            for elmt in array[element+1:]:
                if elmt == array[element]: count += 1
            counts.append(count)

    return unique_array, counts