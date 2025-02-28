from pathlib import Path
import netCDF4                  # can use this or scipy module in line below
import numpy as np
import matplotlib.pyplot as plt

datafolderpath = Path('hyades_output/109103/')
laserData = np.loadtxt(datafolderpath/'laser_pulse_shape_final.txt', delimiter=',') # laser pulse shape
datafile = netCDF4.Dataset(datafolderpath/'hyChDD.cdf') # .cdf output file

# Importing Laser Pulse Shape
laserTime = laserData[:,0]
laserPow  = laserData[:,1]

# Variables
    # R = zone boundary positions (cm) 
    # tn = Fusion energy emission rate
    # fE = Fusion energy released
    # Tr = Radiation temperatures
    # dene = Electron density
variables_list = ['R','Rcm','Rho','Ti','Te','Pres','Bpeprdr','Bpeprd','Tr','Dene','DumpTimes']
variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']

data = dict((label, datafile.variables[var][:]) for label, var in zip(variable_labels, variables_list))

# Plotting R vs Time
xdata = 1E9*data['time'] # times (now nanoseconds),
ydata = 1E4*data['r'] # Zone boundaries (um). Transposed such that each index is one zone boundary from 0 to 1058.

fig, ax = plt.subplots(figsize=(8,5)) 

# Plotting (a sample of 100) zone boundaries 
indices = np.linspace(0, len(xdata)-1, 200).astype(int)
print(indices)
ax.plot(xdata[indices], ydata[indices], c="black", lw=0.2)
ax.plot([], [], c="black", lw=0.2, label='Zone Boundary') # Dummy Plot for Legend

ax2 = ax.twinx()
ax2.plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, zorder=10, label="Laser Pulse")

ax.set_ylim([0,500])
ax.set_ylabel('Radius ($\\mathrm{\\mu}$m)')
ax.set_xlabel('Time (ns)')
plt.show()