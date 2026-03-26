from pathlib import Path
import netCDF4                  
import numpy as np
import matplotlib.pyplot as plt
import sys

from src.main import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, RadialProfile, PrintTimes, RadialProfileSlider

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

    # === Line Plots ===
    # LinePlot_Radius_v_Time(data, laserTime, laserPow, highlight_zones=[79,227])  # Plot zone boundary positions vs time, with highlighted zones
    # for time in [2.65,2.75,2.85,2.95,3.05]:    # Plot MMI acquisition times on top of plotted graph
    #     plt.axvline(x=time, color='blue', linestyle='--', lw=1, alpha=0.5)  

    # # Snapshot of radial profile at 2.9ns
    # time = 3
    # fig, ax = plt.subplots(figsize=(6,3))
    # ax2 = ax.twinx()
    # Te_line = RadialProfile(data, time=time, ax=ax, color="#3072b1", title=None, xlim=(-120,120), plot_shell_boundary=False)
    # ne_line = RadialProfile(data, time=time, variable="dene", ylabel="$n_e$ ($10^{{24}}$ cm$^{-3}$)", color="#A72626", label="$n_e$", linestyle="--", title=f"$T_e$ and $n_e$ Radial Profile t = {time} ns", xlim=(-120,120), ax=ax2)
    # lines = Te_line + ne_line
    # labels = [l.get_label() for l in lines]
    # plt.legend(lines, labels, loc=0)
    # # plt.show()

    # Radial profile of Te and ne with slider to select time
    fig, ax_te, ax_dene, slider = RadialProfileSlider(data, time=2.9, xlim=(-200,200))

    ## === Color Plots (x,y,z = time,radius,`variable`) === 
    # Colormap(data, laserTime, laserPow, variable="te", log=False)
    Colormap(data, laserTime, laserPow, variable="rho", log=True, highlight_zones=[78,227])

    plt.show()
