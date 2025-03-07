compile_opt idl2

filepath= 
ncid=ncdf_open("./hyades_output/109103/hyChDD.cdf",/nowrite)

;get a variable out of the data
NCDF_VARGET, ncid,  varid, var_name

;Close the netcdf file
NCDF_CLOSE, ncid

end