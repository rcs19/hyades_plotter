#!/bin/bash
#SBATCH --job-name=dipB_delay0,56dipB_dt0,08
#SBATCH --output=hyades_output.txt
#SBATCH -p demagnete
#SBATCH --ntasks=1
#SBATCH --time=04:59:59
srun hyades /work3/clf/rscott/hyades/implosions/direct/omega/big_dip/210113a_CH14DT45Rmp20Pk28_1,9d0/dipB_delay0,56dipB_dt0,08/hyImp35.inf
