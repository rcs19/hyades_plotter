from matplotlib import pyplot as plt
import matplotlib

size=12
matplotlib.rc('font', size=size)
matplotlib.rc('axes', titlesize=size)

xdata = [1.7, 1.8, 1.9]
ydata = [1.5, 2.0, 0.8]

fig, ax = plt.subplots()    
ax.scatter(xdata, ydata, s=100)
ax.set_ylabel("DD Neutron Yield ($\\times 10^{10}$)")
ax.set_xlabel("Shock Timing (ns)")
ax.grid(alpha=0.2)
ax.set_title("Neutron Yields")
plt.show()