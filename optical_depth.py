import numpy as np
import matplotlib 
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from pathlib import Path
from src.plotfuncs import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, RadialProfile, PrintTimes, RadialProfileSlider, GetArealDensity, GetMass, GetWeightedAvg, PlotArealDensity, PlotAvgVar
from scipy.interpolate import RegularGridInterpolator

def read_int_opacity(filepath):
    """
    Reads a .prp file and extracts temperature, density, combined mesh,
    and integrated Rosseland mean opacity values from the
    "Int. Rosseland Mean Opacity  (cm**2/g)" section.

    Returns:
        temp_grid (np.ndarray): 1D array of length 22 (Temperatures)
        density_grid (np.ndarray): 1D array of length 20 (Densities)
        combined_grid (np.ndarray): 2D array of shape (22, 20, 2) containing (Temp, Density) pairs
        rosseland_opacity (np.ndarray): 2D array of shape (22, 20) containing opacity values
    """
    header_target = "Int. Rosseland Mean Opacity"
    
    temp_grid = []
    density_grid = []
    opacity_flat = []
    
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # 1. Parse Mesh Parameters for Opacity (nT=22, nD=20)
    for i, line in enumerate(lines):
        if "mesh  parameters for opacity" in line:
            # Read nT and temperature grid
            nT = int(lines[i + 4].strip())

            # Temperatures span across multiple lines until nD is reached
            idx = i + 5
            temp_vals = []
            while len(temp_vals) < nT:
                temp_vals.extend([float(x) for x in lines[idx].split()])
                idx += 1
            temp_grid = np.array(temp_vals)

            # Read nD and density grid
            nD = int(lines[idx].strip())
            idx += 1
            den_vals = []
            while len(den_vals) < nD:
                den_vals.extend([float(x) for x in lines[idx].split()])
                idx += 1
            density_grid = np.array(den_vals)
            break

    # 2. Parse Integrated Rosseland Mean Opacity block (22 x 20 = 440 values)
    target_found = False
    for i, line in enumerate(lines):
        if header_target in line:
            target_found = True
            idx = i + 1
            total_elements = len(temp_grid) * len(density_grid)
            while len(opacity_flat) < total_elements and idx < len(lines):
                curr_line = lines[idx].strip()
                # Stop if a new header section begins
                if curr_line.startswith("*****") or "Planck Mean Opacity" in curr_line:
                    break
                opacity_flat.extend([float(x) for x in curr_line.split()])
                idx += 1
            break

    if not target_found:
        raise ValueError(f"Header '{header_target}' not found in the file.")

    # Now we have an array of 440 opacity values. This is [n1t1, n1t2, ..., n1t22, n2t1, n2t2, ..., n20t22] in row-major order. Reshape into 20 density rows x 22 temperature cols
    rosseland_opacity = np.array(opacity_flat).reshape(len(density_grid),len(temp_grid))

    # Construct the combined 20 rows x 22 columns mesh grid containing (Temperature, Density) tuples
    D_mesh, T_mesh = np.meshgrid(density_grid, temp_grid) #indexing='ij')
    combined_grid = np.stack((D_mesh, T_mesh), axis=-1)  # Shape: (20, 22, 2)
    
    return np.array(density_grid), np.array(temp_grid), combined_grid, rosseland_opacity

def read_opacity(filepath):
    """
    Returns a 3D array containing k(nu) for each temperature and density.
    [[n1, T1], [n2, T2], ... [n20, T22]] - each one of these arrays [n1, T1] is a 1D array of length 500 for the 500 photon energy groups.
    """
    temp_grid       = []
    density_grid    = []
    group_structure = []

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        # 1. Parse Mesh Parameters for Opacity (nT=22, nD=20)
        if "mesh  parameters for opacity" in line:
            # Read nT and temperature grid
            nT = int(lines[i + 4].strip())

            # Temperatures span across multiple lines until nD is reached
            idx = i + 5
            temp_vals = []
            while len(temp_vals) < nT:
                # this reads the temperature block until the number of temperature values is reached
                temp_vals.extend([float(x) for x in lines[idx].split()])
                idx += 1
            temp_grid = np.array(temp_vals)

            # Read nD and density grid
            nD = int(lines[idx].strip())
            idx += 1
            den_vals = []
            # this reads the density block until the number of density values is reached
            while len(den_vals) < nD:
                den_vals.extend([float(x) for x in lines[idx].split()])
                idx += 1
            density_grid = np.array(den_vals)

            # after density grid is number of photon energy groups
            n_group = int(lines[idx].strip())

        #2. Parse Photon Group Structure for Opacity 
        if "group structure (eV)" in line:
            idx = i + 1
            group_vals = []
            while len(group_vals) < n_group+1: # there should be n_groups+1 values for group structure as these define lower and upper bounds of each group
                group_vals.extend([float(x) for x in lines[idx].split()])
                idx += 1
            group_structure = np.array(group_vals)
            break

    grid_total = nT * nD
    if grid_total != (len(temp_grid)*len(density_grid)):
        raise ValueError(f"Mismatch in expected grid sizes: nT*nD={nT*nD} vs len(temp_grid)*len(density_grid)={len(temp_grid)*len(density_grid)}")
    
    # now there should be grid_total + 1 blocks with header containing "Rosseland Mean Opacity (cm**2/g)". The first one is "Int. Rosseland Mean Opacity (cm**2/g)" but we can skip this.
    # compile all 440 blocks of 500 values into a 3D array of shape (nD, nT, n_group) = (20, 22, 500)
    opacity_table = []
    for i, line in enumerate(lines):
        if "Rosseland Mean Opacity (cm**2/g)" in line and "Int." not in line:
            idx = i + 1
            kappa_nu = []
            while len(kappa_nu) < n_group and idx < len(lines):
                curr_line = lines[idx].strip()
                # Stop if a new header section begins - but the block should break after n_group values are read
                if curr_line.startswith("*****") or "Planck Mean Opacity" in curr_line:
                    print(f"Warning: Reached unexpected header before reading all n_group={n_group} opacity values at line {i}.")
                    break
                kappa_nu.extend([float(x) for x in curr_line.split()])
                idx += 1
            if len(kappa_nu) != n_group:
                raise ValueError(f"Expected {n_group} opacity values, but got {len(kappa_nu)} at block starting at line {i}.")
            opacity_table.append(kappa_nu)
            # group_structure = np.vstack((group_structure, opacity_flat)) if group_structure.size else np.array(opacity_flat)
    opacity_table = np.array(opacity_table)
    
    # len(k_nu) = 500, and the flat opacity_table will be [k_nu(n1, T1), k_nu(n1, T2), ..., k_nu(n20, T22)] of shape (nD*nT, n_group). 
    # Reshape opacity_table into (nD, nT, n_group) = (20, 22, 500)
    opacity_table = opacity_table.reshape(len(density_grid), len(temp_grid), n_group)  # Shape: (nD, nT, n_group)
    print(f"opacity_table shape: {opacity_table.shape}") # should be (nD*nT, n_group) = (440, 500)
    print(f"Done loading opacity data: {len(density_grid)} densities, {len(temp_grid)} temperatures, {n_group} photon energy groups.")

    D_mesh, T_mesh = np.meshgrid(density_grid, temp_grid) #indexing='ij')
    combined_grid = np.stack((D_mesh, T_mesh), axis=-1)  # Shape: (20, 22, 2)

    return np.array(density_grid), np.array(temp_grid), combined_grid, group_structure, opacity_table

def interpolate_opacity(n, T, filename=None, temp_grid=None, density_grid=None, opacity_grid=None):
    """
    Interpolates integrated Rosseland mean opacity for a given Temperature (T)
    and Density/Number Density (n) using log-log bilinear interpolation.

    Parameters:
        T (float): Temperature value
        n (float): Density / Number density value
        filename (str): Path to the opacity data file
        temp_grid (np.ndarray): Optional pre-loaded temperature grid
        density_grid (np.ndarray): Optional pre-loaded density grid
        opacity_grid (np.ndarray): Optional pre-loaded opacity grid

    Returns:
        float: Interpolated opacity value (linear scale)
    """
    # Load grids and opacity table
    if temp_grid is None or density_grid is None or opacity_grid is None:
        temp_grid, density_grid, _, opacity_grid = read_int_opacity(filename)

    # Validate inputs against domain bounds
    if not (temp_grid.min() <= T <= temp_grid.max()):
        raise ValueError(f"Temperature T={T} out of grid bounds [{temp_grid.min()}, {temp_grid.max()}]")
    if not (density_grid.min() <= n <= density_grid.max()):
        raise ValueError(f"Density n={n} out of grid bounds [{density_grid.min()}, {density_grid.max()}]")

    # Transform coordinates and values to log10 space
    log_T_grid = np.log10(temp_grid)
    log_n_grid = np.log10(density_grid)
    
    # Avoid log of zero or negative opacity values if any exist
    log_opacity_grid = np.log10(np.maximum(opacity_grid, 1e-30))

    # Construct 2D bilinear interpolator on log grid
    interpolator = RegularGridInterpolator(
        (log_n_grid, log_T_grid),
        log_opacity_grid,
        method="linear",
        bounds_error=True
    )

    # Evaluate point in log space
    point_log = np.array([np.log10(n), np.log10(T)])
    log_opacity_val = interpolator(point_log)[0]

    # Convert back to linear scale
    return 10**log_opacity_val

def plot_opacity_grid(filename="opacity_data.txt", num_points=200):
    """
    Plots a 2D color map of the integrated Rosseland opacity grid with 
    density (x-axis) and temperature (y-axis) on log scales using 
    bilinear interpolation.

    Parameters:
        filename (str): Path to the opacity data file
        num_points (int): Resolution of the grid for smooth bilinear rendering
    """
    # 1. Load data
    temp_grid, density_grid, combined_grid, opacity_grid = read_int_opacity(filename)

    # Define the interpolator once and evaluate on the mesh instead of `interpolate_opacity()` which would repeatedly create the interpolator
    log_T_grid = np.log10(temp_grid)
    log_n_grid = np.log10(density_grid)

    interpolator = RegularGridInterpolator(
    (log_n_grid, log_T_grid),
    opacity_grid,
    method="linear",
    bounds_error=True
    )

    tmin = 1e1
    tmax = 2000
    nmin = 1e20
    nmax = 5e25

    # 2. Create fine logarithmic evaluation grid
    log_T_fine = np.linspace(np.log10(tmin), np.log10(tmax), num_points)
    log_n_fine = np.linspace(np.log10(nmin), np.log10(nmax), num_points)
    T_fine = 10**log_T_fine
    n_fine = 10**log_n_fine

    X_n, Y_T = np.meshgrid(log_n_fine, log_T_fine)

    # Evaluate the interpolator on the fine mesh
    opacity_fine = interpolator((X_n, Y_T))

    # 4. Generate Plot
    plt.figure(figsize=(9, 7))
    
    # Render fine mesh with pcolormesh (shading='bilinear' or 'gouraud')
    pcm = plt.pcolormesh(
        n_fine, 
        T_fine, 
        opacity_fine, 
        norm=LogNorm(vmin=opacity_fine.min(), vmax=opacity_fine.max()),
        cmap='inferno',
        shading='gouraud'
    )

    # Overlay discrete original data points
    orig_N, orig_T = np.meshgrid(density_grid, temp_grid)
    plt.scatter(orig_N, orig_T, color='cyan', s=10, alpha=0.5, label='Grid Points')

    # Formatting axes and labels
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Density $n$ (${g/cm}^3$)')
    plt.ylabel('Temperature $T$ (K)')
    plt.title('Integrated Rosseland Mean Opacity Map')
    plt.xlim(nmin, nmax)
    plt.ylim(tmin, tmax)
    
    # Add Colorbar
    cbar = plt.colorbar(pcm)
    cbar.set_label('Opacity $\\kappa$ ($\\text{cm}^2/\\text{g}$)')
    
    plt.legend(loc='upper right')
    plt.tight_layout()

if __name__=='__main__':
    density_grid, temp_grid, combined_grid, group_structure, opacity = read_opacity("CH_opacity_study.prp")
    # # get index of photon energy closest to 3000 eV in group_structure
    # idx_i = np.argmin(np.abs(group_structure - 3000))
    # idx_f = np.argmin(np.abs(group_structure - 5000))

    # # Option 1: Sum opacity from 3000 to 5000, or idx_i to idx_f, for each density and temperature pair. This will give us a 2D array of shape (nD, nT) = (20, 22)
    # # opacity = np.sum(opacity[:,:,idx_i:idx_f], axis=2) # sum over the last axis (photon energy groups)

    # # Option 2: Only take opacity at a specific photon energy
    # opacity = opacity[:,:,idx_f] 

    if True:
        # plot_opacity_grid(filename="CH_opacity_study.prp", num_points=300)
        # density_grid, temp_grid, combined_grid, opacity = read_int_opacity("CH_opacity_study.prp")

        ## === Load in Data ===
        datafolderpath = Path('shots/98263/')
        data, laserTime, laserPow = Load_Data(datafolderpath)
        time = 2.9484 # 2.9484

        # extract temperature and densities at this time
        # get time index of closest time
        # print(data["te"].shape) # this should be something like (time, cells) = (468, 377)
        time_index = np.argmin(np.abs(data["time"] - time*1e-9))
        r_400_idx = np.argmin(np.abs(data["r"][time_index,:]*1e4 - 400))
        print(f"ne at 400um: {data['dene'][time_index, r_400_idx]:.2e} cm^-3, Te at 400um: {data['te'][time_index, r_400_idx]*1e3:.2f} eV")

        # slice spatial indices for relevant shell region only! 
        min_index = 79
        max_index = r_400_idx
        dr = np.diff(data["r"][time_index,:])
        rhoR     = data["rho"][time_index, min_index:] * dr[min_index:]
        r_slice  = data["rcm"] [time_index, min_index:max_index]*1e4
        te_slice = data["te"]  [time_index, min_index:max_index]*1e3 # keV to eV
        ne_slice = data["dene"][time_index, min_index:max_index]

        # interpolate opacity for each specified photon energy
        log_ne = np.log10(ne_slice)
        log_te = np.log10(te_slice)
        log_ne_te = list(zip(log_ne, log_te))         # combine these into (ne, te) for each cell for the interpolator e.g. [(n1,t1), (n2,t2), ... (n377,t377)]

        fig, ax = plt.subplots(nrows=3, figsize=(6,6), sharex=True)#, layout="tight")
        ax_te   = ax[0]
        ax_rhoR = ax[0].twinx()
        ax_op   = ax[1]
        ax_tau  = ax[2]
        ax0_lines = []
        fig.subplots_adjust(hspace=0.1)

        # draw vertical lines at shell boundaries for each subplot
        for axes in ax:
            axes.axvline(x=data["r"][time_index, 78]*1e4, color='black', linestyle='--', lw=1, alpha=0.5)
            axes.axvline(x=data["r"][time_index, 227]*1e4, color='black', linestyle='--', lw=1, alpha=0.5)

        ax_ne = ax_te.twinx()
        Te_line = RadialProfile(data, time=time, ax=ax_te, color="#CC2D2D", label="$\\mathrm{T_e}$", xlabel=None, title=None, xlim=(-140,140), plot_shell_boundary=False)
        ne_line = RadialProfile(data, time=time, variable="dene", ylabel="$\\mathrm{n_e}$ ($\\times 10^{{24}}$ cm$^{-3}$)", color="#009721", label="$\\mathrm{n_e}$", linestyle="--", title=f"", xlim=(-120,120), ax=ax_ne, plot_shell_boundary=False) # title="$T_e$ and $n_e$ Radial Profile, t = {time} ns"
        # put ne on the left y-axis
        ax_te.set_yticks(np.arange(0, 2.1, 0.5))
        ax_ne.yaxis.set_label_position("left")
        ax_ne.yaxis.set_ticks_position("left")
        ax_ne.spines['left'].set_position(('outward', 60))

        lines = Te_line + ne_line

        ax_tau.set_ylabel('Optical Depth $\\tau$')

        ax_op.set_ylabel('Opacity $\\mathrm{\\kappa}$ ($\\mathrm{cm^2/g}$)')
        # ax_op.spines['right'].set_position(('outward', 60))
        ax_op.tick_params(axis='y')

        energies = [5000, 10000]
        cmap = plt.get_cmap('jet')
        # use colormap to get colors for each energy

        for i, hnu in enumerate(energies):
            # first = 0, last = 1
            color = cmap(i/(len(energies)-1)) if len(energies) > 1 else cmap(0.5)
            # color = cmap(i/len(energies)) # evenly spaced colors for each energy
            op_hnu = opacity[:,:,np.argmin(np.abs(group_structure - hnu))] # get opacity at this photon energy
            opacity_interpolator = RegularGridInterpolator((np.log10(density_grid), np.log10(temp_grid)), np.log10(op_hnu), method='linear', bounds_error=True)
            log_opacity = opacity_interpolator(log_ne_te) # now for each cell, we can compute the opacity
            tau = 10**log_opacity * rhoR[:max_index-min_index]
            # plot optical depth
            ln_tau = ax_tau.plot(r_slice, np.cumsum(tau), color=color, label=f'${hnu*1e-3:.0f}$ keV', linewidth=2) 
            # plot opacity
            op_line = ax_op.plot(r_slice, 10**log_opacity, color=color, label=f'${hnu*1e-3:.0f}$ keV', linewidth=2) 
            
        ax_rhoR.set_ylabel('Shell $\\mathrm{\\rho R}$ (mg/cm$^2$)')
        ln_rhoR = ax_rhoR.plot(data["r"][time_index,min_index:-1]*1e4, np.cumsum(rhoR)*1e3, label='$\\mathrm{\\rho R}$', linestyle='-')
        lines += ln_rhoR

        labels = [l.get_label() for l in lines]
        ax0_labels = [l.get_label() for l in ax0_lines]
        ax[0].legend(lines, labels, loc="upper right")
        ax[1].legend(loc="upper left")
        ax[2].legend(loc="upper left")
        ax[2].set_xlim(0,150)    
        ax[2].set_xlabel("Radius (µm)")

        plt.show()