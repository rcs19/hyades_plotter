from scipy.io import netcdf_file
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

filepath = Path('hyades_output/SAI/hyImp35.cdf')
data = netcdf_file(filepath) # Load in the output from the simulation
for var in data.variables.keys():
    print(var)