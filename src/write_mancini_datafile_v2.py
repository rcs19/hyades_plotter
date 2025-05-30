"""
Authors: Ryan Saputil, Robbie Scott (original script)
---
Script which converts the .cdf output file from a Hyades simulation into a .txt file readable by the post-processor code written by Dr. Roberto Mancini (UNR).
This verison no longer relies on the `pyh2d` package and the header of the output file should now be in the correct format.  


How to use:
---
1. cd {directory containing .cdf output file}  
2. python3 write_mancini_datafile_v2.py
This should output a .txt file in the correct format named `rmancini_hyades_data_{current directory name}`

Alternatively you can call the function `write_mancini_file` with the path to the .cdf file as an argument. This is useful for batch converting items. This would look something like:
```
for shotfolder in shotdir:
    write_mancini_file(shotfolder / 'hyChDD.cdf')
```
    

Notes:
---
1. The header title and region labels are currently hardcoded in...
2. The header of the file should look like this:

Mancini hyades:860OD, 27um CH, fill 20atm D, 0.25%Ar
  1 D 0.25%Ar          2   79
  2 CH                80  227
  3 CH               227  378

The last two numbers of each line represent the first and last Lagrangian zones of that region (inclusive both ends)! The first zone of the first region is always `2` to make it compatible with the post-processor.


TO-DO:
---
Code in the total neutrons:
    # total neutrons = dndt[current time] * timestep for all times  
    nts = 0.0
    dt = times[1]
        for ts in range(ntimes):
            nts += (n_rate_t[ts]*1e-9*dt)

"""

import os
import datetime
import numpy as np
from netCDF4 import Dataset
from pathlib import Path

def write_mancini_file(cdf_file):
    """
    Write a file readable by Dr. Roberto Mancini's post-processor. 
    Arguments:
        cdf_file (pathlike): path to the .cdf hyades output file.
    
    """
    # --- Load Data ---
    with Dataset(cdf_file, mode='r') as f:
        reg = f.variables['RegNums'][1:-1]
        atm_num = f.variables['AtmNum'][:]
        atm_wgt = f.variables['AtmWgt'][:]
        atm_frc = f.variables['AtmFrc'][:]
        times = f.variables['DumpTimes'][:] * 1e9
        r = f.variables['R'][:] * 1e+4
        vr = f.variables['U'][:] * 1e-5
        den = f.variables['Rho'][:]
        tele = f.variables['Te'][:]
        tion = f.variables['Ti'][:]
        ne = f.variables['Dene'][:]
        bp_pwr = f.variables['Bpeprdr'][:] * 1e-7
        bp_eng = f.variables['Bpeprd'][:] * 1e-7
        bp_edep = f.variables['Bpedep'][:] * 1e-7
        vol = f.variables['Vol'][:, 1:]
        mass = vol * den

        try:
            ne_iso = f.variables['Tndeni'][:]
            ne_iso_t = ne_iso.T
        except KeyError:
            print('TN isotope densities not found')
            ne_iso_t = None

        try:
            ubin0 = f.variables['Ubin'][:] * 1e-7
            grp_bnds = f.variables['PhotonGroupBnds'][:]
            grp_avgs = f.variables['PhotonGroupAvgs'][:]
            ubin = 0.5 * (ubin0[:, :-1] + ubin0[:, 1:])  # Assuming center averaging
            ubin = ubin * vol[:, :, np.newaxis]
        except KeyError:
            print('Variable ubin failed to read in')
            ubin = None
            grp_bnds = None
            grp_avgs = None

    # --- Initialize Arrays ---
    nr = reg.shape[0]
    nz = 2
    ntimes = times.size

    reg_arr = np.zeros((nr, nz), dtype=int)
    r_t = np.zeros((nr + 1, nz, ntimes))
    vr_t = np.zeros((nr + 1, nz, ntimes))
    den_t = np.zeros((nr, nz, ntimes))
    tele_t = np.zeros((nr, nz, ntimes))
    tion_t = np.zeros((nr, nz, ntimes))
    bp_pwr_t = np.zeros((nr, nz, ntimes))
    ne_t = np.zeros((nr, nz, ntimes))

    # --- Populate Arrays ---
    reg_arr[:, 1] = reg
    r_t[:, 1, :] = r.T
    vr_t[:, 1, :] = vr.T
    den_t[:, 1, :] = den.T
    tele_t[:, 1, :] = tele.T
    tion_t[:, 1, :] = tion.T
    bp_pwr_t[:, 1, :] = bp_pwr.T

    try:
        ne_t[:, 1, :] = ne.T
    except ValueError:
        print("Failed to read in electron density: RT calculations won't work!")

    # --- Neutron Production Rate ---
    n_rate_t = 0.8 * np.sum(bp_pwr_t[:, 1, :], axis=0) / (14.1e6 * 1.6e-19)  # n/sec

    # --- Output File Creation ---
    shotnum = os.path.basename(os.getcwd())
    # now = datetime.datetime.now()
    # yymmdd = f"{now.year % 100:02d}{now.month:02d}{now.day:02d}"
    # print(yymmdd)

    output_filename = f'rmancini_hyades_data_{shotnum}.txt'

    reg_names = ['D 0.25%Ar', 'CH', 'CH']
    num_regions = np.max(reg_arr)

    with open(output_filename, 'w') as file:
        file.write('Mancini hyades:860OD, 27um CH, fill 20atm D, 0.25%Ar\n')

        for reg_num in range(num_regions):
            reg_cells = np.where(reg_arr[:, 1] == reg_num + 1)[0]
            if reg_cells.size > 0:
                min_cell = np.min(reg_cells)+2 # Note: first zone of first region must be 2 (needed for post-processor)
                max_cell = np.max(reg_cells)+2
                file.write(f'{reg_num + 1:3d} {reg_names[reg_num]:<15}{min_cell:5d}{max_cell:5d}\n')  
                # Reg name is 16 characters including space behind it4
                # reg_cells is zone number (inclusive)
                # for example region 1: 1 to 78, region 2: 79 to 226 inclusive on both sides.
    # Original  file.write('%3.0d %s %3.0d %3.0d\n'%(reg+1,reg_names[reg],np.min(reg_cells),np.max(reg_cells)))        

        file.write('Distance (cm),Elec. density (cm^-3),Elec. temp. (eV), Mass density (g/cm3), Velocity (cm/s), Ion tem. (eV)\n')

        for ts in range(ntimes):
            file.write(f'time= {times[ts]:2.4f}     dndt(1/s)= {n_rate_t[ts]:2.4e}     nts= 0.000E+00\n')
    #       file.write('time= %2.4f     dndt(1/s)= %2.4e     nts= 0.000E+00\n'%(h2d.g.times[ts],h2d.g.n_rate_t[ts]))    # Original
            for ir in range(r_t.shape[0]):
                try:
                    data_line = (
                        r_t[ir, 1, ts] * 1e-4,
                        ne_t[ir, 1, ts],
                        tele_t[ir, 1, ts] * 1e3,
                        den_t[ir, 1, ts],
                        vr_t[ir, 1, ts] * 1e5,
                        tion_t[ir, 1, ts] * 1e3
                    )
                except IndexError:
                    data_line = (
                        r_t[ir, 1, ts] * 1e-4,
                        0, 0, 0,
                        vr_t[ir, 1, ts] * 1e5,
                        0
                    )
                file.write('{:15.5e} {:15.5e} {:15.5e} {:15.5e} {:15.5e} {:15.5e}\n'.format(*data_line))

    print(f'\n(ntimes, nzones)={ne.shape}')
    print(f"Wrote file '{output_filename}'")

if __name__=="__main__":
    """
    If this script is executed in a shot foler, it will look for a .cdf file `hyChDD.cdf` and convert this for the post-processor.
    """
    current_dir = Path.cwd()
    cdf_file = current_dir / 'hyChDD.cdf'
    write_mancini_file(cdf_file)