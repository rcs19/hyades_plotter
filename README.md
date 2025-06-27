# Hyades_Plotter

This repository has code for plotting NetCDF output files from Hyades. The plots can either be created in Python/Jupyter (`.py` and `.ipynb` files) or in IDL (`.pro` files).

## Python Script - `main.py`

This script contains plotting functions for Hyades outputs and requires a .cdf NetCDF output file. If you do not have a .cdf file, you can use the Hyades utility `PPF2NCDF` to generate one from the .ppf (see Hyades manual).

## IDL Script - `main_idl.pro`

Starting IDL from command line on YPI servers (I am using Sausage):
```
$ /opt/york/phys/pkg/Modules/idl/idl71/bin/idl
```
To run the IDL script:
```
IDL> .run ./src/main_idl.pro
```

### Notes:

- In the Hyades input deck `hyChDD.inf` the laser intensity is in ergs cm$^{-2}$ s$^{-1}$. This corresponds to $10^{-7}$ W cm$^{-2}$ (1 erg = $10^{-7}$ J). In `laser_pulse_shape_final.txt` the units of time are ns and laser intensity is $10^{12}$ cm$^{-2}$.


