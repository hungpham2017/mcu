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


def make_basis_dict(cell):
    '''Give a PySCF cell object, return the a dict of basis fucnctions'''
    ao_labels = cell.ao_labels()
    species = []
    lm_list = []
    for ao in ao_labels:
        atm, orb = ao.split()[1:]
        species.append(atm)
        lm_list.append(orb.strip()[1:])
    return species, lm_list
 
