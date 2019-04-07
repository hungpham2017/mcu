# Modeling and Crystallographic Utilities (mcu)
- A package for post periodic wave function and crystallography analysis.
- **mcu** is designed for large scale analysis rather than a code for generating figures for one calculation.
- Currently **mcu** only supports VASP outputs but other electronic codes are under development and will be released soon.

<img src="https://github.com/hungpham2017/mcu/blob/master/doc/MoS2.png" width="500" align="middle">
Projected band structure of MoS2 using the color map style

# Prerequisites
- The code was tested for the below libraries, the older versions, however, can work too. Just need to be tested. I do it soon ...
- Python 3.5 or higher
- Numby 1.15.4
- Matlibplot 3.0.1

 
# Current functions:
- A class to extract all the info from vasprun.xml 
- Functions to read WAVEDER (binary) and WAVEDERF (formatted)
- Plotting band structure, projected band structure with different styles (will be documented soon)
- Computing bandgap


# Future functions:
 - Optical adsorption coefficient
 - Set up a EOS calculation
 - Phonon calculation
 - Spin texture
 - Band structure unfolding
 - Topological analysis
 - Connecting to Bilbao Crystallographic Server
