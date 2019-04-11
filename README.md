# Modeling and Crystallographic Utilities (mcu)
A package for post periodic wave function and crystallography analysis. **mcu** is designed for large scale analysis rather than for generating figures for one calculation.

## Getting Started
A laziest demonstration where most of default values are used

```
import mcu
run = mcu.VASP()        # Define a VASP object
mcu.get_bandgap()       # Compute bandgap
mcu.plot_pband()        # Plot projected band structure
```

<img src="https://github.com/hungpham2017/mcu/blob/master/doc/MoS2.png" width="500" align="middle">
Projected band structure of MoS2 using the color map style

## Prerequisites
- The code was tested for the below libraries, the older versions, however, can work too. Just need to be tested. I do it soon ...
- Python 3.5 or higher
- Numby 1.15.4
- Matlibplot 3.0.1

 
## Current functions:
- Currently **mcu** only supports VASP outputs but other electronic codes are under development and will be released soon.
- A class to extract all the info from vasprun.xml. All the infomation (starting/intermediate/final structure, electronic inputs,...)
can be easily accessed via the **vasprun** attribute of an **mcu.VASP** object. Will be documented soon.
- Functions to read WAVEDER (binary) and WAVEDERF (formatted)
- Plotting band structure, projected band structure with different styles (will be documented soon)
- Computing bandgap

## Function will be released soon:
 - Spin texture
 
## Future functions:
 - Optical adsorption coefficient
 - Set up a EOS calculation
 - Phonon calculation
 - Band structure unfolding
 - Topological analysis
 - Connecting to Bilbao Crystallographic Server
 
## Bugs and suggesting useful functions:
There is a function that you would like it to be included in **mcu**. Shoot me an email or open an issue here!
I don't promise to put everything you like in **mcu** because sometimes you can modify the code with ease for your need.
Necessary functions for sure will be included.


## Authors
- **Hung Pham** - [My Github](https://github.com/hungpham2017)

## License
This project is licensed under the Apache License - see the [LICENSE](LICENSE) file for details