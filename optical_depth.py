import numpy as np
import matplotlib 
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm
from pathlib import Path
from src.plotfuncs import Load_Data, LinePlot_Radius_v_Time, LinePlot_v_Time, Colormap, RadialProfile, PrintTimes, RadialProfileSlider, GetArealDensity, GetMass, GetWeightedAvg, PlotArealDensity, PlotAvgVar
from scipy.interpolate import RegularGridInterpolator

def read_opacity_table(filename):
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
    
    with open(filename, 'r') as f:
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

    # Reshape the flat opacity list into a 22x20 array
    rosseland_opacity = np.array(opacity_flat).reshape(len(temp_grid), len(density_grid))

    # Construct the combined 22x20 mesh grid containing (Temperature, Density) tuples
    T_mesh, D_mesh = np.meshgrid(temp_grid, density_grid, indexing='ij')
    combined_grid = np.stack((T_mesh, D_mesh), axis=-1)

    return np.array(temp_grid), np.array(density_grid), combined_grid, rosseland_opacity

def interpolate_opacity(T, n, filename="opacity_data.txt"):
    """
    Interpolates integrated Rosseland mean opacity for a given Temperature (T)
    and Density/Number Density (n) using log-log bilinear interpolation.

    Parameters:
        T (float): Temperature value
        n (float): Density / Number density value
        filename (str): Path to the opacity data file

    Returns:
        float: Interpolated opacity value (linear scale)
    """
    # Load grids and opacity table
    temp_grid, density_grid, _, opacity_grid = read_opacity_table(filename)

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
        (log_T_grid, log_n_grid),
        log_opacity_grid,
        method="linear",
        bounds_error=True
    )

    # Evaluate point in log space
    point_log = np.array([np.log10(T), np.log10(n)])
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
    temp_grid, density_grid, _, opacity_grid = read_opacity_table(filename)

    # 2. Create fine logarithmic evaluation grid
    log_T_fine = np.linspace(np.log10(temp_grid.min()), np.log10(temp_grid.max()), num_points)
    log_n_fine = np.linspace(np.log10(density_grid.min()), np.log10(density_grid.max()), num_points)
    
    n_mesh, T_mesh = np.meshgrid(10**log_n_fine, 10**log_T_fine)
    opacity_fine = np.zeros_like(T_mesh)

    # 3. Evaluate interpolated values across the 2D mesh
    for i in range(num_points):
        for j in range(num_points):
            opacity_fine[i, j] = interpolate_opacity(T_mesh[i, j], n_mesh[i, j], filename=filename)

    # 4. Generate Plot
    plt.figure(figsize=(9, 7))
    
    # Render fine mesh with pcolormesh (shading='bilinear' or 'gouraud')
    pcm = plt.pcolormesh(
        n_mesh, 
        T_mesh, 
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
    plt.xlabel(r'Density $n$ ($\text{g/cm}^3$)', fontsize=12)
    plt.ylabel(r'Temperature $T$ (K)', fontsize=12)
    plt.title('Integrated Rosseland Mean Opacity Map', fontsize=14)
    
    # Add Colorbar
    cbar = plt.colorbar(pcm)
    cbar.set_label(r'Opacity $\kappa$ ($\text{cm}^2/\text{g}$)', fontsize=12)
    
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    # --- Example Usage ---
    # temp_grid, density_grid, combined_grid, opacity = read_opacity_table("CH_finer.prp")
    # print("Temperature Grid shape:", temp_grid.shape)        # (22,)
    # print("Density Grid shape:", density_grid.shape)            # (20,)
    # print("Combined Grid shape:", combined_grid.shape)          # (22, 20, 2)
    # print("Rosseland Opacity Array shape:", opacity.shape)      # (22, 20)    

    kappa = interpolate_opacity(T=20, n=6e25, filename="CH_finer.prp")
    print("Interpolated Opacity:", kappa)

    plot_opacity_grid(filename="CH_finer.prp", num_points=500)