import numpy as np
from matplotlib import pyplot as plt
plt.rcParams.update({'font.size': 12})
from pathlib import Path

def load_spec(filepath):
    spec = np.loadtxt(filepath, skiprows=4, comments='#')
    spec = spec.T
    return spec

coreonly_pth = Path("spect3d\\20260113_Y2TAP1_98263_2.9ns\\98263_f3_coreonly_spec.dat")
coreshell_pth = Path("spect3d\\20260113_Y2TAP1_98263_2.9ns\\98263_f3_coreshell_spec.dat")

coreonly = load_spec(coreonly_pth)
coreshell = load_spec(coreshell_pth)

divide = coreshell[1]/coreonly[1]

fig, ax = plt.subplots()
ln_core = ax.plot(coreonly[0], coreonly[1], label='Core Only', color="#3e6bce")
ln_coreshell = ax.plot(coreshell[0], coreshell[1], label='Core + Shell', color="#77C965")
ax.set_xlabel('Energy (eV)')
ax.set_ylabel("Detector Flux (erg/cm$^{2}$/eV)")

ax2 = ax.twinx()
ax2.set_ylabel("$\\approx$ Shell Transmission Factor")
ln_divide = ax2.plot(coreonly[0], divide, label='(Core + Shell) ÷ (Core Only)', color="#a02c2c", linestyle='--')

lns = ln_core + ln_coreshell + ln_divide
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='center right')
# ax.set_yscale('log')
ax.set_title("Spect3D 98263 2.9ns - Core Only vs Core + Shell")
plt.show()