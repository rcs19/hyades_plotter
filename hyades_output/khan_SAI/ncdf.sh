#!/bin/bash
#SBATCH --job-name=ncdf_convesion
#SBATCH --output=ppf2ncdf.txt
#SBATCH -p demagnete
#SBATCH --ntasks=1
#SBATCH --time=05:00
ppf2ncdf hyImp35.ppf 
