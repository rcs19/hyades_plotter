from pathlib import Path
import netCDF4                  
import numpy as np
import matplotlib.pyplot as plt
import sys

from src.main import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, RadialProfile, PrintTimes

if __name__=='__main__':

    """
    For reference: variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']
    All variables use CGS units except temperatures which are in keV, or as otherwise noted.

    See Hyades User Manual p. 45 for variables and units
    - fE (BPEPRD)   : fusion Energy, zone TN burn energy produced, erg
    - tn (BPEPRDR)  : zone thermonuclear burn energy production rate, erg/sec
    - tr (TR)       : zone radiation temperature, keV
    - dene (DENE)   : zone electron density, cm^-3  
    """
    ## === Load in Data ===
    datafolderpath = Path('shots/98263/')
    data, laserTime, laserPow = Load_Data(datafolderpath)

    # # === Line Plots ===
    LinePlot_Radius_v_Time(data, laserTime, laserPow) 
    for time in [2.65,2.75,2.85,2.95,3.05]:    # Plot MMI acquisition times on top of plotted graph
        plt.axvline(x=time, color='blue', linestyle='--', lw=1, alpha=0.5)  

    # LinePlot_v_Time(data, laserTime, laserPow, variable='te') 
    # fig, ax = plt.subplots()
    # l1 = RadialProfile(data, time=3.0, variable="dene", ylabel="$n_e$ (cm$^{-3}$)", title="Electron Density Radial Profile", xlim=(-100,100), ax=ax)
    # # ax2 = ax.twinx()
    # # l2 = RadialProfile(data, time=3.0, variable="rho", ylabel="$\\rho$ (g cm$^{-3}$)", title="", xlim=(-100,100), ax=ax, color='orange')
    # ax3 = ax.twinx()
    # # l3 = RadialProfile(data, time=3.0, variable="te", ylabel="$T_e$ (keV)", title="", xlim=(-100,100), ax=ax2, color='green')
    # l4 = RadialProfile(data, time=3.0, variable="p", ylabel="p (Gbar)", title="", xlim=(-100,100), ax=ax3, color='blue')
    # lns = l1 + l4
    # labs = [l.get_label() for l in lns]
    # ax.legend(lns, labs, loc=0)

    ## === Color Plots (x,y,z = time,radius,`variable`) === 
    # Colormap(data, laserTime, laserPow, variable="te", log=False)
    Colormap(data, laserTime, laserPow, variable="te", log=False)

    plt.show()
