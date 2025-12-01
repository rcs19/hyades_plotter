from pathlib import Path
import netCDF4                  # can use this or scipy module in line below
import numpy as np
import matplotlib.pyplot as plt
import sys

from src.main import Load_Data, RadialProfile

def load_dti(filepath):
    """ Load Spect3D .dti output file and extract Radius and Intensity columns"""

    # Extract Radius (2nd column) and Intensity (3rd column)
    data = np.loadtxt(filepath, comments='#')
    radius = data[:, 1]
    intensity = data[:, 2]

    xdata = np.concatenate([-radius[::-1], radius])*1e4  # Convert from cm to microns
    ydata = np.concatenate([intensity[::-1], intensity])

    image = np.array([xdata, ydata])

    return image

if __name__=='__main__':

    """
    For reference: variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']
    See Hyades User Manual p. 45 for variables and units
    """
    ## === Load in Data ===
    datafolderpath = Path('shots/98263/')
    data, laserTime, laserPow = Load_Data(datafolderpath)

    times = data['time'] * 1E9 # Convert to ns

    # === Radial Profile at given time ===
    time = 2.85 # Edit this (ns)
    xlim = (-120, 120)  # Edit this (microns)

    for time in [2.751,2.85,2.95,3.05]:
        tindex = np.argmin(np.abs(times - time))  # Find timedump index at this time 
        print(f"{tindex}, {times[tindex]:.4g}")

        dti_path = Path(f"spect3d/20251201_98263_fuelonly_t268-345/results_01/t0{tindex}/20251201_98263_fuelonly_t268-345_0{tindex}.dti")

        fig, ax = plt.subplots(nrows=2, figsize=(7,4), sharex=True)

        # Top Plot: Radial Profile of Intensity
        intensity = load_dti(dti_path)
        ax[0].plot(intensity[0], np.log10(intensity[1]), color="#00b45a")
        ax[0].set_ylabel("$log_{{10}}$(I) (arb.)")
        ax[0].set_title(f"(Top) Spect3D Intensity Profile\n(Bottom) $T_e$ and $n_e$ Profile")
        ax[0].text(0.95, 0.95, f'$t = {time:.2f}$ ns',horizontalalignment='right',verticalalignment='top',transform = ax[0].transAxes)
        # Bottom Plot: Radial Profile of Te and ne
        ax2 = ax[1].twinx()
        Te_line = RadialProfile(data, time=time, color="#3072b1", title=None, xlim=xlim, ax=ax[1])
        ne_line = RadialProfile(data, time=time, variable="dene", ylabel="$n_e$ ($\\times 10^{24}$ cm$^{-3}$)", color="#A72626", label="$n_e$", linestyle="--", title="", xlim=xlim, ax=ax2)
        lines = Te_line + ne_line
        labels = [l.get_label() for l in lines]
        plt.legend(lines, labels, loc=0)
        # fig.savefig('hyades_radial_profile.svg', format="svg", bbox_inches="tight")
        
        fig.subplots_adjust(hspace=0)

    plt.show()
