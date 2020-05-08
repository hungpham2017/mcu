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


# source: https://www.science.co.il/elements/
ELEMENTS = {
'H':[1,1.008],
'D':[1,2.014],
'He':[2,4.003],
'Li':[3,6.941],
'Be':[4,9.012],
'B':[5,10.811],
'C':[6,12.011],
'N':[7,14.007],
'O':[8,15.999],
'F':[9,18.998],
'Ne':[10,20.180],
'Na':[11,22.990],
'Mg':[12,24.305],
'Al':[13,26.982],
'Si':[14,28.086],
'P':[15,30.974],
'S':[16,32.065],
'Cl':[17,35.453],
'Ar':[18,39.948],
'K':[19,39.098],
'Ca':[20,40.078],
'Sc':[21,44.956],
'Ti':[22,47.867],
'V':[23,50.942],
'Cr':[24,51.996],
'Mn':[25,54.938],
'Fe':[26,55.845],
'Co':[27,58.933],
'Ni':[28,58.693],
'Cu':[29,63.546],
'Zn':[30,65.390],
'Ga':[31,69.723],
'Ge':[32,72.640],
'As':[33,74.922],
'Se':[34,78.960],
'Br':[35,79.904],
'Kr':[36,83.800],
'Rb':[37,85.468],
'Sr':[38,87.620],
'Y':[39,88.906],
'Zr':[40,91.224],
'Nb':[41,92.906],
'Mo':[42,95.940],
'Tc':[43,98.000],
'Ru':[44,101.070],
'Rh':[45,102.906],
'Pd':[46,106.420],
'Ag':[47,107.868],
'Cd':[48,112.411],
'In':[49,114.818],
'Sn':[50,118.710],
'Sb':[51,121.760],
'Te':[52,127.600],
'I':[53,126.905],
'Xe':[54,131.293],
'Cs':[55,132.906],
'Ba':[56,137.327],
'La':[57,138.906],
'Ce':[58,140.116],
'Pr':[59,140.908],
'Nd':[60,144.240],
'Pm':[61,145.000],
'Sm':[62,150.360],
'Eu':[63,151.964],
'Gd':[64,157.250],
'Tb':[65,158.925],
'Dy':[66,162.500],
'Ho':[67,164.930],
'Er':[68,167.259],
'Tm':[69,168.934],
'Yb':[70,173.040],
'Lu':[71,174.967],
'Hf':[72,178.490],
'Ta':[73,180.948],
'W':[74,183.840],
'Re':[75,186.207],
'Os':[76,190.230],
'Ir':[77,192.217],
'Pt':[78,195.078],
'Au':[79,196.967],
'Hg':[80,200.590],
'Tl':[81,204.383],
'Pb':[82,207.200],
'Bi':[83,208.980],
'Po':[84,209.000],
'At':[85,210.000],
'Rn':[86,222.000],
'Fr':[87,223.000],
'Ra':[88,226.000],
'Ac':[89,227.000],
'Th':[90,232.038],
'Pa':[91,231.036],
'U':[92,238.029],
'Np':[93,237.000],
'Pu':[94,244.000],
'Am':[95,243.000],
'Cm':[96,247.000],
'Bk':[97,247.000],
'Cf':[98,251.000],
'Es':[99,252.000],
'Fm':[100,257.000],
'Md':[101,258.000],
'No':[102,259.000],
'Lr':[103,262.000],
'Rf':[104,261.000],
'Db':[105,262.000],
'Sg':[106,266.000],
'Bh':[107,264.000],
'Hs':[108,277.000],
'Mt':[109,268.000]
}

spacegroup =[]

