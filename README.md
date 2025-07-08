# Hyades_Plotter

This repository has code for plotting NetCDF output files from Hyades. The plots can either be created in Python/Jupyter (`.py` and `.ipynb` files) or in IDL (`.pro` files).

## Python Script - `main.py`

This script contains plotting functions for Hyades outputs and requires a .cdf NetCDF output file. If you do not have a .cdf file, you can use the Hyades utility `PPF2NCDF` to generate one from the .ppf (see Hyades manual).

This script can:
1. Plot the time evolution of a specified variable for each Lagrangian zone.
2. Plot the time evolution of a specified variable in a colourmap, with radius on the y-axis and time on the x-axis and the variable on the z-(colour)-axis.
3. Plot a radial profile of a specified variable at a specified time.

### Notes:

- In the Hyades input deck `hyChDD.inf` the laser intensity is in ergs cm$^{-2}$ s$^{-1}$. This corresponds to $10^{-7}$ W cm$^{-2}$ (1 erg = $10^{-7}$ J). In `laser_pulse_shape_final.txt` the units of time are ns and laser intensity is $10^{12}$ cm$^{-2}$.


