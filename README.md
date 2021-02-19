# MCMC_atmo
Emcee based code to study planetary atmospheres

The code creates a posterior class, that contains a model, a prior and a likelihood.

Each iteration draws a prior following the prior probability, creates a model and evaluates its likelihood.

A launch script for titan is provided in launch.sh
