import numpy as np
import matplotlib 
import matplotlib.pyplot as plt
from pathlib import Path

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
    data_a, laserTime, laserPow = Load_Data(Path('shots/109092'))
    data_b, laserTime, laserPow = Load_Data(Path('shots/109092Ti'))

    # === Line Plots ===
    # LinePlot_Radius_v_Time(data_a, laserTime, laserPow, highlight_zones=[79,227])  # Plot zone boundary positions vs time, with highlighted zones
    # for time in [2.65,2.75,2.85,2.95,3.05]:    # Plot MMI acquisition times on top of plotted graph
    #     plt.axvline(x=time, color='blue', linestyle='--', lw=1, alpha=0.5)  


    # Radial profile of Te and ne with slider to select time
    # time = 2.851
    # time = 2.907
    # time = 3.0
    time = 2.987

    # time = 2.964

    # fig, ax_te, ax_dene, slider = RadialProfileSlider(data_a, time=time, xlim=(-150,150), ymax=[2.6, 16])
    
    if False:
        fig, ax = plt.subplots(figsize=(5,3))
        ax2 = ax.twinx()
        Te_line = RadialProfile(data_a, time=time, ax=ax, color="#CC2D2D", title=None, xlim=(-140,140), plot_shell_boundary=False)
        # rho_line = RadialProfile(data_a, time=time, variable="p", ylabel="$\\rho$ (g/cm$^3$)", color="#A72626", label="$\\rho$", linestyle="--", title=f"$T_e$ and $\\rho$ Radial Profile, t = {time} ns", xlim=(-120,120), ax=ax2)
        ne_line = RadialProfile(data_a, time=time, variable="dene", ylabel="$n_e$ ($\\times 10^{{24}}$ cm$^{-3}$)", color="#009721", label="$n_e$", linestyle="--", title=f"$T_e$ and $n_e$ Radial Profile, t = {time} ns", xlim=(-120,120), ax=ax2, plot_shell_boundary=False)
        # lines = Te_line + rho_line
        lines = Te_line + ne_line
        labels = [l.get_label() for l in lines]
        ax.set_xlim(-2,125)    
        plt.legend(lines, labels, loc="upper right")
        ax.set_xlabel("Radius (µm)")

    # ## === Color Plots (x,y,z = time,radius,`variable`) === 
    # # Colormap(data_a, laserTime, laserPow, variable="te", log=False)
    # font size
    matplotlib.rcParams['font.size'] = 12
    axcolormap_a    = Colormap(data_a, laserTime, laserPow, variable="dene", log=True, highlight_zones=[78])
    axcolormap_b = Colormap(data_b, laserTime, laserPow, variable="dene", log=True, highlight_zones=[78])
    axcolormap_a.set_xlim(2.5,3.5)
    axcolormap_b.set_xlim(2.5,3.5)
    axcolormap_a.set_ylim(0,90)
    axcolormap_b.set_ylim(0,90)

    
    # axcolormap   .set_title("CH DD")
    # axcolormap_Ar.set_title("CH DD Ar")
    # GetArealDensity(data_a, time=time, r1=78, r2=377)
    # # GetArealDensity(data_a, time=time, r1=78, r2=227)
    # GetMass(data_a, time=time, r1=78, r2=130)
    # GetWeightedAvg(data_a, var="te", weightvar="mass", time=time, r1=0, r2=78)
    # GetWeightedAvg(data_b, var="te", weightvar="mass", time=time, r1=0, r2=78)
    # GetWeightedAvg(data_a, var="dene", weightvar="rhoR", time=time, r1=78, r2=130)
    # GetWeightedAvg(data_a, var="dene", weightvar="rhoR", time=time, r1=130, r2=227)
    # Colormap(data_a, laserTime, laserPow, variable="p", log=False, highlight_zones=[78])

    print(np.sum(data_a["fE"][-1,:]))
    print(np.sum(data_b["fE"][-1,:]))
    print(np.sum(data_a["fE"][-1,:])/np.sum(data_b["fE"][-1,:]))

    print()
    plt.show()

