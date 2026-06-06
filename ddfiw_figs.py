from pathlib import Path
import pandas as pd
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

    # time = 2.93

    # fig, ax_te, ax_dene, slider = RadialProfileSlider(data, time=time, xlim=(0,100), ymax=[2.3, 1e25])
    

    downscatter_rhoR = np.array([[2.95, 145, 6], [2.95, 118, 10]]) 
    summary = pd.read_csv("98252.csv", header=0, sep=",")
    summary_times_and_r = pd.read_csv("98252_TIMING.csv", header=0, sep=",")
    # summary_times_and_r["est_time"] = summary_times_and_r["pulse_time"] + 0.06

    # join by both "tim" and "frame" columns
    summary = pd.merge(summary, summary_times_and_r, on=["tim", "frame"], how="inner")
    # print(summary[["tim","frame","est_time","r","r_err"]].to_csv(index=False))

    # SHIFTS and UNITS
    # =======
    data["time"] = data["time"] - 50e-12 # bang time shift
    summary["est_time"] = summary["est_time"] + 0.06
    # summary["ne_avg"] = summary["ne_avg"] * 1e-24
    # summary["ne_std"] = summary["ne_std"] * 1e-24

    fig, ax = plt.subplots(1, figsize=(4, 4))

    for tshift, tim, marker in zip([-2950,-2950,-2950], [3,4,5], ["^", "o", "X"]):
        df_tim = summary[summary["tim"] == tim]
        this_time = -(df_tim["est_time"]*1e3 + tshift)
        ax.errorbar(this_time, df_tim["rhor_avg"]*1e3, yerr=df_tim["rhor_std"]*1e3, fmt=marker, label=f"MMI TIM{tim}", 
                    markersize=10, capsize=0,  
                    mfc="#00FFF2", mec="#006071", ecolor="#006071") # , xerr=np.ones_like(df_tim["rhor_avg"].values)*50
        
    ax.errorbar(-(downscatter_rhoR[0, 0]*1e3 - 2950), downscatter_rhoR[0, 1], yerr=downscatter_rhoR[0, 2], fmt="s", label="WRF1", 
                markersize=10, capsize=0,
                mfc="#FFA8D9", mec="#71003E", ecolor="#71003E")
    ax.errorbar(-(downscatter_rhoR[1, 0]*1e3 - 2950), downscatter_rhoR[1, 1], yerr=downscatter_rhoR[1, 2], fmt="v", label="WRF2", 
                markersize=10, capsize=0,
                mfc="#FFA8D9", mec="#71003E", ecolor="#71003E")
    ax.legend() 
    ax.set_title("Shell Areal Density vs Time")
    ax.set_xlabel("$t-t_{bang}$ (ps)")
    ax.set_xlim(-30,330)
    ax.set_ylabel("$\\mathrm{\\rho R}$ (mg/cm$^2$)")
    ax.xaxis.set_inverted(True)

    # ax.set_xticks([0, 100, 200,300])
    plt.tight_layout()
    # plt.show()

    kbf_data = pd.read_csv("98252_kbframed.txt", comment="#", skiprows=1, header=0, sep="\\s+")
    r_time = data["time"]*1e9 # time array across all zones
    r78 = data["r"][:,78]*1e4 # radius of zone 78 across all times

    # Plotting
    # ========
    markersize = 7

    fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(9,10), sharex=True)
    ax[0].errorbar(kbf_data["t"]*1e-3 + 2.65, kbf_data["r0"],yerr=kbf_data["r0_err"], label="KBFramed", color="#CC2D2D", fmt="x", linestyle="", markersize=markersize, mew=2, capsize=0)
    ax[0].plot(r_time,r78, label="Simulation\nCore radius", color="#000000", ls="-")
    ax[0].set_ylabel("HS Radius (µm)")
    ax_dene = ax[1].twinx()
    line_te = PlotAvgVar(data, variable="te", r1=0, r2=20, ax=ax[1], label="$T_e$", color="#CC2D2D", linestyle="-")
    # data["dene"] = data["dene"]*1e-24
    ln_dene = PlotAvgVar(data, variable="dene", r1=0, r2=77, ax=ax_dene, label="$n_e$", color="#009721", linestyle="--")
    ax[1].set_ylabel("$\\mathrm{T_e}$ (keV)")
    ax_dene.set_ylabel("$\\mathrm{n_e}$ ($\\times 10^{24}$ cm$^{-3}$)")
    ax_tn = ax[2].twinx()
    PlotArealDensity(data, r1=78, r2=227, ax=ax[2], label="$\\rho R$")
    PlotAvgVar(data, variable="tn", r1=0, r2=77, ax=ax_tn, label="Simulation\nneutron prod. rate", color="#000000", alpha=0.7, linestyle="-.")

    # Plot MMI data
    for tshift, tim, marker in zip([0, 0.0, 0.0], [3, 4, 5], ["^", "o", "X"]):
        df_tim = summary[summary["tim"] == tim]

        # 1/e Radius
        ax[0].errorbar(df_tim["est_time"]+tshift, df_tim["r"], yerr=df_tim["r_err"], fmt=marker, label=f"MMI{tim}", 
                    markersize=markersize, capsize=0,
                    mfc="#B6B7CF", mec="#37315A", ecolor="#37315A")
        
        # Electron temperature
        ax[1].errorbar(df_tim["est_time"]+tshift, df_tim["te_avg"]*1e-3, yerr=df_tim["te_std"]*1e-3, fmt=marker, label=f"MMI{tim}", 
                    markersize=markersize, capsize=0,  
                    mfc="#FFAE00", mec="#3D2100", ecolor="#3D2100") # xerr=np.ones_like(df_tim["ne_avg"].values)*0.05
        
        # Electron density
        ax_dene.errorbar(df_tim["est_time"]+tshift, df_tim["ne_avg"], yerr=df_tim["ne_std"], fmt=marker, label=f"MMI{tim}", 
                    markersize=markersize, capsize=0,  
                    mfc="#61FF4C", mec="#006928", ecolor="#006928") # , xerr=np.ones_like(df_tim["ne_avg"].values)*0.05

        # Areal density
        ax[2].errorbar(df_tim["est_time"]+tshift, df_tim["rhor_avg"]*1e3, yerr=df_tim["rhor_std"]*1e3, fmt=marker, label=f"MMI{tim}", 
                    markersize=markersize, capsize=0,  
                    mfc="#00FFF2", mec="#006071", ecolor="#006071")
        

    # Plot nTOF and WRF downscatter data
    ax[2].errorbar(downscatter_rhoR[0, 0], downscatter_rhoR[0, 1], yerr=downscatter_rhoR[0, 2], fmt="v", label="WRF1", 
                markersize=markersize, capsize=0,
                mfc="#FFA8D9", mec="#71003E", ecolor="#71003E")
    ax[2].errorbar(downscatter_rhoR[1, 0], downscatter_rhoR[1, 1], yerr=downscatter_rhoR[1, 2], fmt="s", label="WRF2", 
                markersize=10, capsize=0,
                mfc="#FFA8D9", mec="#71003E", ecolor="#71003E")
    
    for i in range(len(ax)):
        # ax[i].axvline(x=2.914, color='gray', linestyle='--', lw=2, alpha=0.5)  # shock convergence
        ax[i].axvline(x=2.95, color='gray', linestyle='--', lw=2, alpha=0.5)  # bang time

    ax[0].legend(loc="upper left")
    ax[1].legend(loc="upper left")
    ax_dene.legend(loc="lower right")
    ax_tn.set_ylabel("$\\mathrm{t_n}$ (erg s$^{-1}$)")
    ax_tn.legend(loc="lower right")
    ax[2].set_ylabel("$\\mathrm{\\rho R}$ (mg/cm$^2$)")
    ax[2].set_xlabel("Time (ns)")
    ax[2].legend(loc="upper left")

    ax[0].set_ylim(16, 60)
    ax[0].set_xlim(2.55,3.07)
    ax[1].set_ylim(0.0, 2.41)
    ax_dene.set_ylim(0, 8e24)
    ax[2].set_ylim(25,190)
    # ax_tn.set_ylim(0, 11e13)

    ax[0].set_title("Hotspot Radius")
    ax[1].set_title("Electron Temperature and Density")
    ax[2].set_title("Areal Density")
    fig.subplots_adjust(hspace=0.25)

    plt.show()

