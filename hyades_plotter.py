from pathlib import Path
import netCDF4                  
import numpy as np
import matplotlib.pyplot as plt
import sys

from src.plotfuncs import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, RadialProfile, PrintTimes, RadialProfileSlider, GetArealDensity, GetMass, GetWeightedAvg

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
    datafolderpath = Path('shots/98252/')
    data, laserTime, laserPow = Load_Data(datafolderpath)

    # === Line Plots ===
    # LinePlot_Radius_v_Time(data, laserTime, laserPow, highlight_zones=[79,227])  # Plot zone boundary positions vs time, with highlighted zones
    # for time in [2.65,2.75,2.85,2.95,3.05]:    # Plot MMI acquisition times on top of plotted graph
    #     plt.axvline(x=time, color='blue', linestyle='--', lw=1, alpha=0.5)  


    # Radial profile of Te and ne with slider to select time
    time = 2.851
    # fig, ax_te, ax_dene, slider = RadialProfileSlider(data, time=time, xlim=(-150,150), ymax=[1.8, 8])
    
    # Snapshot of radial profile at 2.9ns
    time = 2.9
    fig, ax = plt.subplots(figsize=(5,3))
    ax2 = ax.twinx()
    Te_line = RadialProfile(data, time=time, ax=ax, color="#3072b1", title=None, xlim=(-120,120), plot_shell_boundary=False)
    rho_line = RadialProfile(data, time=time, variable="rho", ylabel="$\\rho$ (g/cm$^3$)", color="#A72626", label="$\\rho$", linestyle="--", title=f"$T_e$ and $\\rho$ Radial Profile, t = {time} ns", xlim=(-120,120), ax=ax2)
    lines = Te_line + rho_line
    labels = [l.get_label() for l in lines]
    ax.set_xlim(0,150)    
    plt.legend(lines, labels, loc="upper right")

    ## === Color Plots (x,y,z = time,radius,`variable`) === 
    # Colormap(data, laserTime, laserPow, variable="te", log=False)
    # Colormap(data, laserTime, laserPow, variable="rho", log=True, highlight_zones=[78,227])
    # GetArealDensity(data, time=time, r1=78, r2=232)
    GetArealDensity(data, time=time, r1=78, r2=227)
    GetMass(data, time=time, r1=78, r2=130)
    GetWeightedAvg(data, var="te", weightvar="rhoR", time=time, r1=78, r2=130)
    GetWeightedAvg(data, var="te", weightvar="rhoR", time=time, r1=130, r2=227)

    plt.show()
