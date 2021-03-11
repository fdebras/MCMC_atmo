# MCMC_atmo
Parallel emcee-based code to study planetary atmospheres. 

Python modules needed: dill, schwimmbad, emcee, json, runpy, numpy, corner,  mpi4py, lmfit. Thus far, the models are based on petitRADTRANS (https://petitradtrans.readthedocs.io/en/latest/index.html)

The code creates a posterior class for every walker of the MCMC algorithm, that contains a model, a prior and a likelihood. 

Each iteration draws a prior following the prior probability and the Markhov chain iteration, evaluates the model from the prior parameters and calculates its likelihood.

A launch script for titan is provided in launch.sh
