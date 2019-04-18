[![Build Status](https://travis-ci.com/hungpham2017/mcu.svg?branch=master)](https://travis-ci.com/hungpham2017/mcu)
# Modeling and Crystallographic Utilities (mcu)
A package for periodic wavefunction and crystallography analysis. **mcu** is designed to support large scale analysis and topological descriptions for periodic wavefunction.

## A quick look
A laziest demo where most of default values are used

```
import mcu
run = mcu.VASP()        # Define a VASP object
run.get_bandgap()       # Compute bandgap
run.plot_pband()        # Plot projected band structure
```

<img src="https://github.com/hungpham2017/mcu/blob/master/docs/image/MoS2.png" width="500" align="middle">
Projected band structure of MoS2 using the color map style

## Documentation:
-  Installation and tutorials can be found [here](https://hungpham2017.github.io/mcu)

## Authors
- **Hung Pham** - [My Github](https://github.com/hungpham2017)

## License
This project is licensed under the Apache License - see the [LICENSE](LICENSE) file for details