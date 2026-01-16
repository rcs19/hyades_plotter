from pathlib import Path
import netCDF4                  
import numpy as np
import matplotlib.pyplot as plt
import sys

def Load_Data(datafolderpath):
    """
    Loads in .cdf data and laser data from a Hyades output folder. The output folder usually contains .inf and .ppf files for running the simulation but this script requires the .cdf NetCDF file. If you do not have a .cdf file, you can use the Hyades utility `PPF2NCDF` to generate one from the .ppf (see Hyades manual).
     
    Returns:
        data - dict, contains variables and their corresponding data array
        laserTime - Laser Time data
        laserPow - Laser Power data
    
    """
    try:
        datafile = netCDF4.Dataset(datafolderpath/'hyChDD.cdf') # .cdf output file
    except:
        print("Error: hyChDD.cdf not found.")
        sys.exit()

    try:
        laserData = np.loadtxt(datafolderpath/'laser_pulse_shape_final.txt', delimiter=',') # laser pulse shape
        laserTime = laserData[:,0]        # ns
        laserPow  = laserData[:,1]        # W cm-2
    except:
        laserTime = laserPow = None

    # Variables
        # R = zone boundary positions (cm) 
        # tn = Fusion energy emission rate
        # fE = Fusion energy released
        # Tr = Radiation temperatures
        # dene = Electron density
    variables_list = ['R','Rcm','Rho','Ti','Te','Pres','Bpeprdr','Bpeprd','Tr','Dene','DumpTimes']
    variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']

    data = dict((label, datafile.variables[var][:]) for label, var in zip(variable_labels, variables_list))
    return data, laserTime, laserPow

def LinePlot_Radius_v_Time(data, laserTime, laserPow):
    """
    Plot Zone Boundary Positions vs Time.
    """
    # Plotting R vs Time
    xdata = 1E9*data['time'] # times (now nanoseconds),
    ydata = 1E4*data['r'] # Zone boundaries (um). Transposed such that each index is one zone boundary from 0 to 1058.

    fig, ax = plt.subplots(figsize=(6,4)) 

    gas_boundary = 79
    ch_boundary = 380

    ax.plot(xdata, ydata[:,0:gas_boundary], c="#004B53", lw=0.2)            # Deuterium
    ax.plot(xdata, ydata[:,gas_boundary:ch_boundary], c="#3B3B3B", lw=0.2)  # Plastic

    line1 = ax.plot([], [], c="#3B3B3B", label='CH Plastic', lw=3) # Dummy Plot for Legend
    line2 = ax.plot([], [], c="#004B53", label='Deuterium', lw=3) # Dummy Plot for Legend
    
    lines = line1 + line2
    # Overlay Laser Pulse Shape
    if (laserPow is not None) and (laserTime is not None):
        ax2 = ax.twinx()
        lns3 = ax2.plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, alpha=0.6, zorder=10, label="Laser Pulse")
        lines = lines + lns3
        ax2.set_ylabel("Laser Intensity ($10^{12}$ W cm$^{-2}$)")
    else:
        print("No laser data")

    # --- Click Event Handler ---
    click_marker, = ax.plot([], [], marker='x', color='red', markersize=10, )
    click_text = ax.text(0.95, 0.6, '', fontsize=12, color='white', va='center', ha='right', bbox=dict(facecolor='black', alpha=0.9, edgecolor='none'), transform=ax.transAxes)
    def on_click(event):
        if event.inaxes is not None:
            # 1. Convert click from display pixels to ax's data coordinates
            # This ignores which axes is 'on top' and uses ax's scale
            inv = ax.transData.inverted()
            x_ax, y_ax = inv.transform((event.x, event.y))

            # 2. Find the index of the nearest Time and Radius using transformed coords
            ix = (np.abs(xdata - x_ax)).argmin()
            iy = (np.abs(ydata[ix, :] - y_ax)).argmin()

            # 3. Extract the values
            # Ensure indices are within bounds
            if 0 <= ix < data['te'].shape[0] and 0 <= iy < data['te'].shape[1]:
                te_val = data['te'][ix, iy]      # Removed .T if using raw data indices
                dene_val = data['dene'][ix, iy]
                # 3. Update the marker position
                click_marker.set_data([x_ax], [y_ax])
                # click_text.set_position((x_ax, y_ax+20))
                click_text.set_text(f'T$_{{e}} = {te_val:.3f}$ keV\nn$_{{e}} = {dene_val:.2e}$ g cm$^{{-3}}$')
                fig.canvas.draw_idle() # Redraws everything efficiently                
                print(f"\nClick event: {x_ax:.2f} ns, {y_ax:.2f} µm ({ix},{iy})")
                print(f"Te = {te_val:.3f} keV, ne = {dene_val:.3e} g/cm³")
            else:
                print(data['te'].shape)
                print(f"\nClick event: {ix}, {iy}")
                print("Click was outside the data range.")

    cid = fig.canvas.mpl_connect('button_press_event', on_click)
    
    labels = [l.get_label() for l in lines]
    ax.legend(lines, labels, loc=0)
    ax.set_ylim([0,500])
    ax.set_ylabel('Radius ($\\mathrm{\\mu}$m)')
    ax.set_xlabel('Time (ns)')
    ax.set_title('Time Plot of Lagrangian Zone Boundaries')
    # fig.savefig('r_vs_t_plot.pdf', format='svg') # Optional save figure as pdf vector graphic for Latex figures

def LinePlot_v_Time(data, laserTime, laserPow, variable="r"):
    """
    Plot Zone Variable vs Time
    """
    # Plotting R vs Time
    xdata = 1E9*data['time'] # times (now nanoseconds),
    ydata = data[variable] 

    fig, ax = plt.subplots(figsize=(8,5)) 

    # Plotting (a sample of n) zone boundaries 
    indices = np.linspace(0, len(xdata)-1, 200).astype(int)
    ax.plot(xdata[indices], ydata[indices], c="black", lw=0.2)
    ax.plot([], [], c="black", lw=0.2, label=variable) # Dummy Plot for Legend

    # Overlay Laser Pulse Shape
    if (laserPow is not None) and (laserTime is not None):
        ax2 = ax.twinx()
        ax2.plot(laserTime, laserPow, color="red", linestyle="--", alpha=0.6, linewidth=2, zorder=10, label="Laser Pulse")
        ax2.legend(loc='lower left')
        ax2.set_ylabel("Laser Power (TW)")
    else:
        print("No laser data")

    # ax.set_ylim([0,500])
    ax.set_ylabel(variable)
    ax.set_xlabel('Time (ns)')

def Colormap(data, laserTime, laserPow, variable="dene", log=False):
    """
    Plots the specified variable as a colormap over time and radius.
    """
    # Create Grid
    xdata = 1E9 * data['time']           # Time in nanoseconds
    # y-axis: Radius in micrometers
    if variable == "tr":    # Zone radiation temperature indexes from 0 to nzones
        ydata = 1E4 * data['rcm'][:,:-1]       
    else:                   # Other zone variables index from 1 to nzones
        ydata = 1E4 * data['rcm'][:,1:-1]       
    X, Y = np.meshgrid(xdata, ydata[0, :])  # Meshgrid for plotting

    if log:
        Z_calc = np.log10(data[variable]).T             # Transpose to match meshgrid orientation
        Z = np.where(data[variable].T <= 0, np.min(Z_calc), Z_calc)  # Need this to handle log(0)
    else:
        Z = data[variable].T

    minimum = np.min(Z[0,:])
    # Create figure and axis
    fig, ax = plt.subplots(figsize=(9,5))

    # Main pcolormesh plot
    pc = ax.pcolormesh(X[:,], 
                       ydata.T[:,], 
                       Z,
                       vmin=minimum, # vmin defines the minimum data range, i.e. minimum density
                       cmap='viridis'
                       )

    # Optional laser plot
    if (laserPow is not None) and (laserTime is not None):
        ax2 = ax.twinx()
        ax2.plot(laserTime, laserPow, color="red", linestyle="--", linewidth=2, zorder=10, label="Laser Pulse", alpha=0.6)
        # ax2.legend(loc="lower left")
        ax2.set_ylabel("Laser Power (TW)")
    else:
        print("No Laser Pulse data passed")

    # --- Click Event Handler ---
    click_marker, = ax.plot([], [], marker='x', color='red', markersize=10, )
    click_text = ax.text(0.95, 1.05, '', fontsize=11, color='white', va='center', ha='right', bbox=dict(facecolor='black', alpha=0.9, edgecolor='none'), transform=ax.transAxes)
    def on_click(event):
        if event.inaxes is not None:
            # 1. Convert click from display pixels to ax's data coordinates
            # This ignores which axes is 'on top' and uses ax's scale
            inv = ax.transData.inverted()
            x_ax, y_ax = inv.transform((event.x, event.y))

            # 2. Find the index of the nearest Time and Radius using transformed coords
            ix = (np.abs(xdata - x_ax)).argmin()
            iy = (np.abs(ydata[ix, :] - y_ax)).argmin()

            # 3. Extract the values
            # Ensure indices are within bounds
            if 0 <= ix < data['te'].shape[0] and 0 <= iy < data['te'].shape[1]:
                te_val = data['te'][ix, iy]      # Removed .T if using raw data indices
                dene_val = data['dene'][ix, iy]
                pres_val = data['p'][ix, iy]
                # 3. Update the marker position
                click_marker.set_data([x_ax], [y_ax])
                # click_text.set_position((x_ax, y_ax+20))
                click_text.set_text(f'({x_ax:.2f} ns, {y_ax:.1f} µm): T$_{{e}} = {te_val:.3f}$ keV, n$_{{e}} = {dene_val:.2e}$ cm$^{{-3}}$, p = {pres_val* 1E-6 * 1E-9:.2f} Gbar')
                fig.canvas.draw_idle() # Redraws everything efficiently                
                print(f"\nClick event: {x_ax:.2f} ns, {y_ax:.2f} µm ({ix},{iy})")
                print(f"({x_ax:.2f} ns, {y_ax:.2f} µm): Te = {te_val:.3f} keV, ne = {dene_val:.3e} cm$^{{-3}}$, p = {pres_val* 1E-6 * 1E-9:.2f} Gbar") # Convert pressure from dyne/cm2 (= 1 barye) to Gbar (1 baryre = 1E-6 bar)
            else:
                print(data['te'].shape)
                print(f"\nClick event: {ix}, {iy}")
                print("Click was outside the data range.")

    cid = fig.canvas.mpl_connect('button_press_event', on_click)
    
    # Colorbar
    plt.subplots_adjust(right=0.8) # Make space for colorbar
    cax = fig.add_axes([ax.get_position().x1+0.08,ax.get_position().y0,0.04,ax.get_position().height])
    plt.colorbar(pc, cax=cax, label=("log$_{10}$ " if log else " ") + variable) # Similar to fig.colorbar(im, cax = cax)

    # Legend and labels
    ax.set_ylim([0, 500])
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r'Radius ($\mu$m)')

    ### End of plotting code

    # ### Print Z value at given coordinates (time,radius) 
    # userinput = input("Enter `x(ns),y(um)` = ").split(",")
    # x, y = float(userinput[0]), float(userinput[1])

    # i = (np.abs(xdata - x)).argmin()
    # j = (np.abs(ydata[0, :] - y)).argmin()
    # print(i,y)
    # try:
    #     value = Z[i, j]
    #     print(f"Z({x:.2f} ns, {y:.2f} µm) = {value:.3g}")
    # except IndexError:
    #     print(f"Coordinate out of range: ({x}, {y})")

def RadialProfile(data, time=3.0, variable="te", title="Electron Temperature Radial Profile", xlabel="x (um)", ylabel="$T_e$ (keV)", xlim=(0,100), label=None, linestyle = "-", color="black", ax = None):
    """
    Plots radial profile of specified variable. By default, plots electron temperature "te". Data is reflected about x=0.
    
    Arguments
    ---
    data: dict, pandas.dataframe
        - Data loaded in from Load_Data() function 
    time: float
        - Time (ns) to plot
    variable: str
        - Variable to plot on y-axis. Default = "te"
    title:
        - Plot title, default = "Electron Temperature Radial Profile"
    xlabel:
        - Plot x-axis label, default = "x (um)"
    ylabel:
        - y-axis label of plot, default = "$T_e$ (keV)"
    xlim:
        - x-axis limits (um), default = (0,100)
    """
    time_index = (np.abs(data['time'] - float(time)*1e-9)).argmin()
    ydata = data[variable][time_index,:]

    # rcm indexes from 0 to nzones+1, first and last index are 0.0
    if variable in  ['rho','ti','te','p','tn','fE','dene',]: # these variables index from 1 to nzones
        xdata = 1e4 * data['rcm'][time_index,1:-1] 
        if variable == "dene":
            ydata = ydata * 1e-24
        elif variable == "p": 
            ydata = ydata * 1E-6 * 1E-9 # Convert pressure from dyne/cm2 (= 1 barye) to Gbar (1 baryre = 1E-6 bar)
            # print("")
    elif variable == "tr":  # Zone radiation temperature indexes from 0 to nzones
        xdata = 1e4 * data['rcm'][time_index,1:-1] 
        ydata = ydata[1:]   # ignore index 0 where tr = 1e-4
    # Mirror x along the y-axis
    xdata = np.concatenate([-xdata[::-1], xdata])
    ydata = np.concatenate([ydata[::-1], ydata])

    ext_axis = True
    if ax is None:
        fig, ax = plt.subplots()
        ext_axis = False
    lineplot = ax.plot(xdata, ydata, label=(variable if label is None else label), ls=linestyle, color=color)
    ax.set_xlim(xlim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return lineplot

def PrintTimes(data):
    """
    Prints dump times in the data.
    """
    times = data['time'] * 1E9 # Convert to ns
    print("Timedumps (from i=1)\nindex,ns:")
    for i, t in enumerate(times):
        print(f"{i+1},{t:.3f}")

if __name__=='__main__':

    """
    For reference: variable_labels = ['r','rcm','rho','ti','te','p','tn','fE','tr','dene','time']
    See Hyades User Manual p. 45 for variables and units
    """
    ## === Load in Data ===
    datafolderpath = Path('shots/98263/')
    data, laserTime, laserPow = Load_Data(datafolderpath)

    xdata = 1E9*data['time'] # times (now nanoseconds),

    # Time snapshot of radial profile of Te and ne
    time=2.95

    fig, ax = plt.subplots(figsize=(6,3))
    ax2 = ax.twinx()
    Te_line = RadialProfile(data, time=time, ax=ax, color="#3072b1", title=None, xlim=(-120,120))
    ne_line = RadialProfile(data, time=time, variable="dene", ylabel="$n_e$ (cm$^{-3}$)", color="#A72626", label="$n_e$", linestyle="--", title=f"$T_e$ and $n_e$ Radial Profile t = {time} ns", xlim=(-120,120), ax=ax2)
    lines = Te_line + ne_line
    labels = [l.get_label() for l in lines]
    plt.legend(lines, labels, loc=0)
    plt.show()
    # # fig.savefig('hyades_radial_profile.svg', format="svg", bbox_inches="tight")

    plt.show()
