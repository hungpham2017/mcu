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


def save_kpts_bands(filename, list_or_tuple_or_filename):
    kpts, bands = list_or_tuple_or_filename
    np.savez(filename, kpts, np.asarray(bands[0]), bands[1])



def load_kpts_bands(filename):
    data = np.load(filename + '.npz')
    return data['arr_0'], (data['arr_1'], data['arr_2'])
    
