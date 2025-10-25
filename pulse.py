from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

files = ["SI_laser_pulse_shape_final.txt", "SAI_laser_pulse_shape_final.txt"]

fig, ax = plt.subplots(figsize=(8,5))

df_si = pd.read_csv(files[0], sep=",", names = ["time","power"])
df_sai = pd.read_csv(files[1], sep=",", names = ["time","power"])

ax.plot(df_si["time"],df_si["power"]*1e12, lw=3, label="SI")
ax.plot(df_sai["time"],df_sai["power"]*1e12, lw=2, label="SAI")

plt.show()