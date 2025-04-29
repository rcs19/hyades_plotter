from pathlib import Path
import netCDF4                  # can use this or scipy module in line below
import numpy as np
import matplotlib.pyplot as plt
import sys

def Load_Data(datafolderpath):
    """
    Loads in .cdf data and laser data from a Hyades output folder. The output folder usually contains .inf and .ppf files for running the simulation but this script requires the .cdf NetCDF file. If you do not have a .cdf file, you can use the Hyades utility `PPF2NCDF` to generate one from the .ppf (see Hyades manual).
     
    Returns:
        data - dict, contains variables and their corresponding data array
        laserTime - Laser Time data
        laserPow - Laser Power data
    
    """
    try:
        datafile = netCDF4.Dataset(datafolderpath/'hyChDD.cdf') # .cdf output file
    except:
        print("Error: hyChDD.cdf not found.")
        sys.exit()

    try:
        laserData = np.loadtxt(datafolderpath/'laser_pulse_shape_final.txt', delimiter=',') # laser pulse shape
        laserTime = laserData[:,0]
        laserPow  = laserData[:,1]
    except:
        laserTime, laserPow = False

    # Variables
        # R = zone boundary positions (cm) 
        # tn = Fusion energy emission rate
        # fE = Fusion energy released
        # Tr = Radiation temperatures
        # dene = Electron density
    variables_list = ['R','Rcm','Rho','Ti','Te','Pres','Bpeprdr','Bpeprd','Tr','Dene','DumpTimes']
    variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']

    data = dict((label, datafile.variables[var][:]) for label, var in zip(variable_labels, variables_list))
    return data, laserTime, laserPow

def Plot_Radius_v_Time(data, laserTime, laserPow):
    """
    Plot Zone Boundary Radius vs Time.
    """
    # Plotting R vs Time
    xdata = 1E9*data['time'] # times (now nanoseconds),
    ydata = 1E4*data['r'] # Zone boundaries (um). Transposed such that each index is one zone boundary from 0 to 1058.

    fig, ax = plt.subplots(figsize=(8,5)) 

    # Plotting (a sample of n) zone boundaries 
    indices = np.linspace(0, len(xdata)-1, 200).astype(int)
    ax.plot(xdata[indices], ydata[indices], c="black", lw=0.2)
    ax.plot([], [], c="black", lw=0.2, label='Zone Boundary') # Dummy Plot for Legend

    # Overlay Laser Pulse Shape
    try:
        ax2 = ax.twinx()
        ax2.plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, zorder=10, label="Laser Pulse")
        ax2.legend(loc='lower left')
        ax2.set_ylabel("Laser Power (TW)")
    except:
        print("No laser data")

    ax.set_ylim([0,500])
    ax.set_ylabel('Radius ($\\mathrm{\\mu}$m)')
    ax.set_xlabel('Time (ns)')
    plt.show()
    
def Plot_Density_pcolormesh(data, laserTime, laserPow):
    """
    Plots the density over time. 
    Note: Radius data has shape (1058, 530) and density data has shape (1057,530) for some weird reason.
    I had to truncate the y-axis of the grid from 1058 down to 1057 because of this. 
    """
    # Create Grid
    xdata = 1E9 * data['time']           # Time in nanoseconds
    ydata = 1E4 * data['r']              # Radius in micrometers
    X, Y = np.meshgrid(xdata, ydata[0, :])  # Meshgrid for plotting
    Z = np.log10(data['rho'].T)          # Transpose to match meshgrid orientation
    minimum = np.min(np.log10(data['rho'][:, 0]))

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
    try:
        ax2 = ax.twinx()
        ax2.plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, zorder=10, label="Laser Pulse")
        ax2.legend(loc="lower left")
        ax2.set_ylabel("Laser Power (TW)")

    except:
        print("No Laser Pulse data passed")

    # Colorbar
    cbar = fig.colorbar(pc, ax=ax)
    cbar.set_label(r"log$_{10}$(Mass Density (g/cm$^3$))")

    # Legend and labels
    ax.set_ylim([0, 500])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r'Radius ($\mu$m)')
    
    # fig.savefig("./si_siImplosion.png", bbox_inches='tight', pad_inches=0.0)
    plt.show()

if __name__=='__main__':
    # Load in Data
    datafolderpath = Path('hyades_output/109103/')
    data, laserTime, laserPow = Load_Data(datafolderpath)
    
    # Plots
    Plot_Radius_v_Time(data, laserTime, laserPow)
    Plot_Density_pcolormesh(data, laserTime, laserPow)