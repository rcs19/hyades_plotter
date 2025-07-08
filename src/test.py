from pathlib import Path
import netCDF4                  # can use this or scipy module in line below
import numpy as np
import matplotlib.pyplot as plt
import sys
datafolderpath = Path('conor/')

try:
    datafile = netCDF4.Dataset(datafolderpath/'hyChDD.cdf') # .cdf output file
except:
    print("Error: hyChDD.cdf not found.")

# Variables
    # R = zone boundary positions (cm) 
    # tn = Fusion energy emission rate
    # fE = Fusion energy released
    # Tr = Radiation temperatures
    # dene = Electron density
variables_list = ['R','Rcm','Rho','Ti','Te','Bpeprdr','Bpeprd','DumpTimes']
variable_labels = ['r','rcm','rho','ti','te','tn','fE','time']

data = dict((label, datafile.variables[var][:]) for label, var in zip(variable_labels, variables_list))

fusionEnergy    = 1e-7*(np.sum(data['fE'][:],axis=(1))[-1]) #J
numFusReactions = fusionEnergy*(1./(17.6e6*1.6e-19))
print(f"{fusionEnergy:.2e} J or {fusionEnergy*1e7:.2e} ergs")
print(f"{numFusReactions:.2e} fusion reactions")

print(f"Fusion energy = {fusionEnergy:.2e}")

def Colormap(data, laserTime, laserPow, variable="rho", log=False):
    """
    Same as above but for a specified variable.
    """
    # Create Grid
    xdata = 1E9 * data['time']           # Time in nanoseconds
    ydata = 1E4 * data['r']              # Radius in micrometers
    X, Y = np.meshgrid(xdata, ydata[0, :])  # Meshgrid for plotting

    if log:
        Z = np.log10(data[variable]).T          # Transpose to match meshgrid orientation
    else:
        Z = data[variable].T

    minimum = np.min(Z[0,:])
    print(np.max(Z))
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(9,5))

    # Main pcolormesh plot
    pc = ax.pcolormesh(X[:-1,], 
                       ydata.T[:-1,], 
                       Z,
                       vmin=minimum, # vmin defines the minimum data range, i.e. minimum density
                       cmap='viridis'
                       )

    # Optional laser plot
    if (laserPow is not None) and (laserTime is not None):
        ax2 = ax.twinx()
        ax2.plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, zorder=10, label="Laser Pulse")
        ax2.legend(loc="lower left")
        ax2.set_ylabel("Laser Power (TW)")
    else:
        print("No Laser Pulse data passed")

    # Colorbar
    plt.subplots_adjust(right=0.8) # Make space for colorbar
    cax = fig.add_axes([ax.get_position().x1+0.08,ax.get_position().y0,0.04,ax.get_position().height])
    plt.colorbar(pc, cax=cax, label=("log$_{10}$ " if log else " ") + variable) # Similar to fig.colorbar(im, cax = cax)

    # Legend and labels
    ax.set_ylim([0, 500])
    ax.set_xlim([0, 100])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r'Radius ($\mu$m)')
    
    # fig.savefig("./si_siImplosion.png", bbox_inches='tight', pad_inches=0.0)

    plt.show()

    ### End of plotting code

    ### Print Z value at given coordinates (time,radius) 
    userinput = input("Enter `x,y` = ").split(",")
    x, y = float(userinput[0]), float(userinput[1])

    i = (np.abs(xdata - x)).argmin()
    j = (np.abs(ydata[0, :] - y)).argmin()
    print(i,y)
    try:
        value = Z[i, j]
        print(f"Z({x:.2f} ns, {y:.2f} µm) = {value:.3g}")
    except IndexError:
        print(f"Coordinate out of range: ({x}, {y})")

Colormap(data, laserTime=None, laserPow=None, )