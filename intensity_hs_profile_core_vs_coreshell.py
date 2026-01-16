from pathlib import Path
import netCDF4                  # can use this or scipy module in line below
import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd

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

def load_dat(filepath):
    # Load the data
    data = pd.read_csv(filepath, sep="\\s+", skiprows=4, skipfooter=1, engine='python', header=None)
    if filepath.name.endswith("_image.dat"):
        data[0] = data[0] * 1e4     # 1. Convert Radius (Column 0) from cm to microns
        mirrored = data.iloc[::-1].copy()
        mirrored[0] = -mirrored[0]
        combined_data = pd.concat([mirrored, data], ignore_index=True)
        
        return combined_data
    else:
        return data

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
    time = 2.9 # Edit this (ns)
    xlim = (-120, 120)  # Edit this (microns)
    
    r_shell = data['r'][306,78]*1e4  # Example: print radius at index 78 in microns

    spec_core_path = Path("spect3d/core_vs_coreshell_instant/98263_f3_coreonly_spec.dat")
    spec_coreshell_path = Path("spect3d/core_vs_coreshell_instant/98263_f3_coreshell_spec.dat")
    image_core_path = Path("spect3d/core_vs_coreshell_instant/98263_f3_coreonly_image.dat")
    image_coreshell_path = Path("spect3d/core_vs_coreshell_instant/98263_f3_coreshell_image.dat")

    spec_core = load_dat(spec_core_path)
    spec_coreshell = load_dat(spec_coreshell_path)
    image_core = load_dat(image_core_path)
    image_coreshell = load_dat(image_coreshell_path)

    crop = (3400,5000)
    spec_core = spec_core[spec_core[0].between(crop[0], crop[1])]
    spec_coreshell = spec_coreshell[spec_coreshell[0].between(crop[0], crop[1])]
    
    fig, ax = plt.subplots(nrows=2, figsize=(7,4), sharex=True)

    # Top Plot: Radial Profile of Intensity
    ax[0].plot(image_core[0], (image_core[1]), color="#00b45a", label="Core Only")
    ax[0].plot(image_coreshell[0], (image_coreshell[1]), color="#b40087", label="Core + Shell")
    ax[0].axvspan(-r_shell, r_shell, color='gray', alpha=0.2, label='Core Region')
    ax[0].legend()
    ax[1].axvspan(-r_shell, r_shell, color='gray', alpha=0.2, label='Core Region')
    ax[0].set_ylabel("Spectral Flux (erg/cm$^2$/eV)")
    ax[0].set_title(f"(Top) Spect3D Intensity Profile\n(Bottom) $T_e$ and $n_e$ Profile")
    ax[0].text(0.95, 0.95, f'$t = {time:.2f}$ ns',horizontalalignment='right',verticalalignment='top',transform = ax[0].transAxes)
    # Bottom Plot: Radial Profile of Te and ne
    ax2 = ax[1].twinx()
    Te_line = RadialProfile(data, time=time, variable="te", color="#3072b1", title=None, xlim=xlim, ax=ax[1])
    # Tr_line = RadialProfile(data, time=time, variable="tr", color="#fffb00", title=None, xlim=xlim, ax=ax[1])
    ne_line = RadialProfile(data, time=time, variable="dene", ylabel="$n_e$ ($\\times 10^{24}$ cm$^{-3}$)", color="#A72626", label="$n_e$", linestyle="--", title="", xlim=xlim, ax=ax2)#
    lines = Te_line + ne_line # + Tr_line
    labels = [l.get_label() for l in lines]
    plt.legend(lines, labels, loc="upper left")
    # fig.savefig('hyades_radial_profile.svg', format="svg", bbox_inches="tight")
    
    fig.subplots_adjust(hspace=0)

    fig, ax = plt.subplots()
    ax2 = ax.twinx()
    coreline = ax.plot(spec_core[0], spec_core[1], color="#00b45a", label="Core Only")
    coreshellline = ax2.plot(spec_coreshell[0], spec_coreshell[1], color="#b40087", label="Core + Shell")
    lns = coreline + coreshellline
    labels = [l.get_label() for l in lns]
    ax.legend(lns, labels)
    ax.set_ylabel("Spectral Flux (erg/cm$^2$/eV)\nCore Only",)   
    ax2.set_ylabel("Core + Shell",)
    ax.set_xlabel("Energy (eV)")
    plt.show()
