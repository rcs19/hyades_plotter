#!/bin/bash
#SBATCH --job-name=109103
#SBATCH --output=hyades_output.txt
#SBATCH -p demagnete
#SBATCH --ntasks=1
#SBATCH --time=23:59:59
srun hyades /work4/clf/rscott/hyades/implosions/dd/omega/om230725/postshot/240513a_flxLmB0,008/109103/hyChDD.inf
