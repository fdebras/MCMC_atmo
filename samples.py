#!/usr/bin/env python3

import numpy as np

from shutil import copyfile

import h5py
import os


def read_safe(fname):
    temp = fname + ".tmp"
    copyfile(fname, temp)

    f = h5py.File(temp, "r")
    n = f["mcmc"].attrs["iteration"]
    s = f["/mcmc/chain"][()]
    lp = f["/mcmc/log_prob"][()]
    f.close()

    os.remove(temp)

    return s[:n, :, :], lp[:n, :]


param_names = [
#        "kappa_IR",
#        "gamma",
#        "T_int",
#        "T_eq",
        "MMR_H2O",
#        "MMR_CO",
#        "MMR_CO2",
        "Kp",
        "Vsys"
]


def pack_theta(d):
    return np.array([d[p] for p in param_names])
