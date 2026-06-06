from pathlib import Path
import netCDF4                  
import numpy as np
import matplotlib.pyplot as plt
import sys

from src.plotfuncs import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, RadialProfile, PrintTimes, RadialProfileSlider, GetArealDensity, GetMass, GetWeightedAvg, PlotArealDensity, PlotAvgVar

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
    # time = 2.851
    time = 2.907
    # time = 2.947
    # time = 2.987

    time = 2.964

    # fig, ax_te, ax_dene, slider = RadialProfileSlider(data, time=time, xlim=(-150,150), ymax=[2.6, 16])
    
    fig, ax = plt.subplots(figsize=(5,3))
    ax2 = ax.twinx()
    Te_line = RadialProfile(data, time=time, ax=ax, color="#CC2D2D", title=None, xlim=(-120,120), plot_shell_boundary=False)
    # rho_line = RadialProfile(data, time=time, variable="rho", ylabel="$\\rho$ (g/cm$^3$)", color="#A72626", label="$\\rho$", linestyle="--", title=f"$T_e$ and $\\rho$ Radial Profile, t = {time} ns", xlim=(-120,120), ax=ax2)
    ne_line = RadialProfile(data, time=time, variable="dene", ylabel="$n_e$ ($\\times 10^{{24}}$ cm$^{-3}$)", color="#009721", label="$n_e$", linestyle="--", title=f"$T_e$ and $n_e$ Radial Profile, t = {time} ns", xlim=(-120,120), ax=ax2, plot_shell_boundary=False)
    # lines = Te_line + rho_line
    lines = Te_line + ne_line
    labels = [l.get_label() for l in lines]
    ax.set_xlim(-80,80)    
    plt.legend(lines, labels, loc="upper right")
    ax.set_xlabel("Radius (µm)")

    # ## === Color Plots (x,y,z = time,radius,`variable`) === 
    # # Colormap(data, laserTime, laserPow, variable="te", log=False)
    import matplotlib 
    # font size
    matplotlib.rcParams['font.size'] = 12
    axcolormap = Colormap(data, laserTime, laserPow, variable="dene", log=True, highlight_zones=[78])
    axcolormap.set_title("1D HYADES")
    # GetArealDensity(data, time=time, r1=78, r2=377)
    # # GetArealDensity(data, time=time, r1=78, r2=227)
    # GetMass(data, time=time, r1=78, r2=130)
    GetWeightedAvg(data, var="te", weightvar="mass", time=time, r1=0, r2=78)
    GetWeightedAvg(data, var="dene", weightvar="rhoR", time=time, r1=78, r2=130)
    GetWeightedAvg(data, var="dene", weightvar="rhoR", time=time, r1=130, r2=227)

    plt.show()
