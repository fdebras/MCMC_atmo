#!/bin/sh
# -N atmo
#SBATCH --nodes=3 --exclusive
#SBATCH --ntasks-per-node=20 
#SBATCH -n 60 
#SBATCH --job-name=mcmc 
#SBATCH -p mem192
#SBATCH -o HD189_comp.out
#SBATCH -e HD189_comp.err
DIR='/home/fdebras/MCMC_atmo/'
source ~/.bashrc
cd $DIR
export OMP_NUM_THREADS=1
mpirun -np  60 python3  mcmc_atm.py --prior priors_HD189.py --like brogi --nwalkers 200 --nsteps 10000 --out HD189_compord.h5 --start-prior --data data_HD189.py
