import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path

datafolderpath = Path('shots/98252/')

laserData = np.loadtxt(datafolderpath/'laser_pulse_shape_initial.txt', delimiter=',') # laser pulse shape
laserTime = laserData[:,0]*1e-3        # ns
laserPow  = laserData[:,1]             # W cm-2

fig, ax = plt.subplots(figsize=(5, 3))

ax.plot(laserTime, laserPow, color="red")

ax.set_xlim(0, 2.7)
ax.set_ylim(0, 0.5)
# ax.set_title("Delivered Laser Pulse")
ax.set_xlabel("Time (ns)")
ax.set_ylabel("Laser Power (TW/beam)")
plt.show()