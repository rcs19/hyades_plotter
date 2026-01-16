import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
from pathlib import Path
from scipy.interpolate import interp1d

def load_spec(filepath):
    spec = np.loadtxt(filepath, skiprows=4, comments='#')
    spec = spec.T
    return spec

coreonly_pth = Path("spect3d\\core_vs_coreshell_instant\\98263_f3_coreonly_spec.dat")
coreshell_pth = Path("spect3d\\core_vs_coreshell_instant\\98263_f3_coreshell_spec.dat")

coreonly = load_spec(coreonly_pth)
coreshell = load_spec(coreshell_pth)

divide = coreshell[1]/coreonly[1]

# Calculating transmission factor

op_table = pd.read_csv("spect3d\\T-631eV_ni-1e24_planckabs.dat", skiprows=4, sep="\\s+", names=['Energy (eV)', 'Opacity cm2/g'], skipfooter=1, engine='python')
op_table ["transmission"] = np.exp(- op_table['Opacity cm2/g'] * 50 * 20*1e-4) # density 10 g/cm3, thickness 20 microns

# df_transmission = pd.read_csv("spect3d\\crxo_CH_20um_15gcc.txt", skiprows=2, sep="\\s+", names=['Energy (eV)', 'transmission'], skipfooter=1, engine='python')
df_transmission = op_table

# Applying transmission to coreonly spectrum
trans_func = interp1d(df_transmission['Energy (eV)'], df_transmission['transmission'], bounds_error=False, fill_value=(1.0, 0.0))
transmission = trans_func(coreonly[0])
attenuated_spectrum = coreonly[1] * transmission

fig, ax = plt.subplots()
# ln_core = ax.plot(coreonly[0], coreonly[1], label='Core Only', color="#3e6bce")
ln_coreshell = ax.plot(coreshell[0], coreshell[1], label='Core + Shell', color="#77C965")
ln_attenuated = ax.plot(coreonly[0], attenuated_spectrum, label='Core Only $\\times$ Transmission', color="#2f4e90", linestyle='--')

ax.set_xlabel('Energy (eV)')
ax.set_ylabel("Detector Flux (erg/cm$^{2}$/eV)")

# ax2 = ax.twinx()
# ax2.set_ylabel("$\\approx$ Shell Transmission Factor")
ax2 = ax.twinx()
# ln_divide = ax2.plot(coreonly[0], divide, label='(Core + Shell) ÷ (Core Only)', color="#a02c2c", linestyle='--')
# ln_trans = ax2.plot(df_transmission['Energy (eV)'], df_transmission['transmission'], label='CRXO Transmission', color="#3f3f3f", linestyle='--')

ax2.set_ylabel("Transmission")

# lns = ln_core + ln_coreshell + ln_divide+ ln_trans + ln_attenuated  
lns = ln_coreshell + ln_attenuated  
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='lower right')
ax.set_xlim(200, 10000)
# ax.set_yscale('log')
ax.set_title("Spect3D 98263 2.9ns - Core Only vs Core + Shell")
plt.show()
