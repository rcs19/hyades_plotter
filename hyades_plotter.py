from pathlib import Path
import netCDF4                  
import numpy as np
import matplotlib.pyplot as plt
import sys

from src.main import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, Colormap_Density, RadialProfile, PrintTimes

if __name__=='__main__':

    """
    For reference: variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']
    See Hyades User Manual p. 45 for variables and units
    """
    ## === Load in Data ===
    datafolderpath = Path('shots/98263/')
    data, laserTime, laserPow = Load_Data(datafolderpath)

    # === Line Plots ===
    LinePlot_Radius_v_Time(data, laserTime, laserPow) 
    for time in [2.65,2.75,2.85,2.95,3.05]:    # Plot MMI acquisition times on top of plotted graph
        plt.axvline(x=time, color='blue', linestyle='--', lw=1, alpha=0.5)  

    # LinePlot_v_Time(data, laserTime, laserPow, variable='te') 
    # RadialProfile(data, time=3.0, variable="dene", ylabel="$n_e$ (cm$^{-3}$)", title="Electron Density Radial Profile", xlim=(-100,100))

    ## === Color Plots (x,y,z = time,radius,`variable`) === 
    # Colormap_Density(data, laserTime, laserPow)
    Colormap(data, laserTime, laserPow, variable="te", log=True)
    # Colormap(data, laserTime, laserPow, variable="dene", log=True)

    plt.show()
