from pathlib import Path
from main import Load_Data
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

def Animate_Zone(data, laserTime, laserPow, fps=20):
    r = data["r"][:, :-1] * 1e4
    Z = np.log10(data["rho"])
    Z_label = "Density"
    vmin = np.min(Z)
    vmax = np.max(Z)
    time_array = data['time'] * 1e9

    # Setting up grid and first frame
    theta = np.linspace(0, 2 * np.pi, 180)
    R_grid, Theta = np.meshgrid(r[0], theta)
    X = R_grid * np.cos(Theta)
    Y = R_grid * np.sin(Theta)
    Z0 = np.tile(Z[0], (len(theta), 1))

    # Initialise Figure
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6, 8), height_ratios=[8, 3])

    # Density plot
    mesh = ax[0].pcolormesh(X, Y, Z0, shading='auto', cmap='magma', vmin=vmin, vmax=vmax-0.5)
    ax[0].set_aspect('equal')
    ax[0].set_xlabel('x (um)')
    ax[0].set_ylabel('y (um)')
    ax[0].set_title('Frame: 0')
    fig.colorbar(mesh, ax=ax[0], label=Z_label)

    # Laser pulse plot
    ax[1].plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, zorder=10, label="Laser Pulse")
    ax[1].set_xlabel("Time (ns)")
    ax[1].set_ylabel("Power")
    ax[1].legend()
    ax[1].set_xlim([np.min(time_array),np.max(time_array)])
    time_marker = ax[1].axvline(0, ls='-', color='b', lw=1, zorder=10)

    for a in ax:
        a.set_anchor('W')

    def update(i):
        # Update density plot
        Zi = np.tile(Z[i], (len(theta), 1))
        mesh.set_array(Zi.ravel())
        ax[0].set_title(f"Frame: {i}  |  Time: {time_array[i]:.3g} ns")

        # Update time marker on laser plot
        time_marker.set_xdata([time_array[i],time_array[i]])
        return mesh, time_marker

    interval = 1000/fps # interval between each frame in miliseconds
    ani = animation.FuncAnimation(fig, update, frames=len(time_array), interval=interval, blit=False)
    # ani.save("density_animation.gif", writer=PillowWriter(fps=fps))

    plt.show()

def Animate_Radius(data, laserTime=None, laserPow=None, fps=20):
    """
    Animate radial zone boundaries as concentric circles over time,
    colored by material (deuterium gas vs CH ablator).
    """
    times = 1e9 * data['time']        # Convert to ns
    radii = 1e4 * data['r']           # Convert to μm
    Nframes, Nzones = radii.shape

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_aspect('equal')
    ax.set_xlim(-500, 500)
    ax.set_ylim(-500, 500)
    ax.set_xlabel('x (μm)')
    ax.set_ylabel('y (μm)')
    ax.set_title(f'Time: {times[0]:.3f} ns')

    # Define zone types
    gas_indices = range(0, 78)        # Zones 1–78
    ablator_indices = range(78, 377)  # Zones 79–377

    # Create colored circles for each zone
    circle_lines = []
    for i in range(Nzones):
        color = "#2e7e6e" if i in gas_indices else "#57358a"
        c = Circle((0, 0), radius=0, fill=False, color=color, lw=0.5)
        circle_lines.append(c)
        ax.add_patch(c)

    # Legend
    legend_elements = [
        Line2D([0], [0], color='#2e7e6e', lw=2, label='Gas Fill'),
        Line2D([0], [0], color='#57358a', lw=2, label='Plastic'),
    ]
    ax.legend(handles=legend_elements, loc='upper left')

    # Optional laser subplot
    if laserTime is not None and laserPow is not None:
        ax_laser = fig.add_axes([0.65, 0.65, 0.3, 0.25])
        ax_laser.plot(laserTime, laserPow, color='red', linestyle='--', linewidth=1.5)
        time_marker, = ax_laser.plot([times[0], times[0]], ax_laser.get_ylim(), color='blue', lw=1)
        ax_laser.set_title("Laser Pulse")
        ax_laser.set_xlabel("Time (ns)")
        ax_laser.set_ylabel("Power")
    else:
        time_marker = None

    def update(frame):
        current_radii = radii[frame]

        for c, r in zip(circle_lines, current_radii):
            c.set_radius(r)

        ax.set_title(f'Time: {times[frame]:.3f} ns')

        if time_marker:
            time_marker.set_xdata([times[frame], times[frame]])
            time_marker.set_ydata(ax_laser.get_ylim())

        return circle_lines + ([time_marker] if time_marker else [])
    
    interval = 1000/fps # interval between each frame in miliseconds
    ani = animation.FuncAnimation(fig, update, frames=Nframes, interval=interval)
    ani.save("radius_animation.gif", writer=PillowWriter(fps=fps))
    # plt.tight_layout()
    # plt.show()

if __name__=='__main__':
    # Load in Data
    datafolderpath = Path('shots/98246/')
    data, laserTime, laserPow = Load_Data(datafolderpath)

    Animate_Zone(data, laserTime, laserPow)
    # Animate_Radius(data, laserTime, laserPow)
