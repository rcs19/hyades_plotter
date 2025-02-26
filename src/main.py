from pathlib import Path
import netCDF4
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

filepath = Path('hyades_output/SAI/hyImp35.cdf')
data = netCDF4.Dataset(filepath)

