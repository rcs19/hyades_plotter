#!/bin/bash
#SBATCH --job-name=dipB_delay0,80
#SBATCH --output=hyades_output.txt
#SBATCH -p demagnete
#SBATCH --ntasks=1
#SBATCH --time=04:59:59
srun hyades /work3/clf/rscott/hyades/implosions/direct/omega/big_dip/210210a_CH14DT45Rmp20Pk220_1,9MJ_retuned/dipB_delay0,80/hyImp35.inf
