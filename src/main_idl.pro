; IDL Script for Plotting Hyades Data
compile_opt idl2

datafolderpath = 'hyades_output/109103/'

; Importing Laser Data
laserFile = datafolderpath + 'laser_pulse_shape_final.txt'
laserData = READ_CSV(laserFile)
laserTime = laserData.field1
laserPow = laserData.field2

; Importing Hyades output file (netCDF4 .cdf file)
; CDF_OPEN - open ncdf file, gives you a file ID handle to interact with
; CDF_VARGET - looks up a variable inside ncdf file and returns its variable ID (`varid = NCDF_VARID(cdfid, varname)`)
; CDF_VARINQ - list variables
; NCDF_VARID - 
cdfFile = datafolderpath + 'hyChDD.cdf'
cdf_id = NCDF_OPEN(cdfFile, /NOWRITE)

; Store each variable from the .cdf as a variable in memory
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'R'), r
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Rcm'), rcm
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Rho'), rho
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Ti'), ti
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Te'), te
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Pres'), p
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Bpeprdr'), tn
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Bpeprd'), fE
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Tr'), tr
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'Dene'), dene
NCDF_VARGET, cdf_id, NCDF_VARID(cdf_id, 'DumpTimes'), time

; Close the .cdf file now that we're done with it
NCDF_CLOSE, cdf_id


; Plotting
xdata = 1E9 * time      ; time in nanoseconds
ydata = 1E4 * r         ; zone boundary positions in microns

; Sample a subset of ydata (e.g. plot 200 out of 1058 lines)
n_samples = 100
n_zones = N_ELEMENTS(ydata[*, 0])
indices = LONG(FINDGEN(n_samples) * (n_zones-1) / (n_samples-1))

; Draw a blank plot to set axis and frame
WINDOW, 0, XSIZE=600, YSIZE=400

PLOT, xdata, ydata[*,0], /NODATA, $
      YRANGE=[0, 500], $
      XTITLE='Time (ns)', $
      YTITLE='Radius (!9m!Xm)', $
      BACKGROUND = 0

; Looping over each sample (one line) 
FOR i = 0, N_ELEMENTS(indices)-1 DO BEGIN
    OPLOT, xdata, ydata[indices[i], *], COLOR=255, THICK=0.1
ENDFOR

END