#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

from shutil import copyfile

import corner
import json
import h5py
import os

import pandas as pd

import samples

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("h5file")
parser.add_argument("--corner")
parser.add_argument("--chains")
parser.add_argument("--logp")
parser.add_argument("--conv")
parser.add_argument("--post-csv")
args = parser.parse_args()

s, lp = samples.read_safe(args.h5file)
n, nwalkers, ndim = s.shape
print(ndim)


prout = np.where(s[:,:,1]<200)
ss = s[prout]

# nskip = 3 * (n // 4)
nskip = 600
if args.corner is not None:
	rng = samples.pack_theta({
	"MMR_H2O": (-8.0, -2.0),
	"Kp": (100.0, 300.0),
	"Vsys": (-20.0, 20.0)})
    # "lambda_max_layered": (0.6, 0.75),
    # "Y_surface": (0.16, 0.22),
    # "Y_intern": (0.2, 0.5),
    # "Y_core": (0.3, 0.7),
    # "Z_surface": (0.005, 0.040),
    # "Z_intern": (0.0, 0.32),
    # "Z_core": (0.1, 0.5),
    # "delta_S_layered": (0.0, 0.007),
    # "omega_slope_l0": (-1, 1),
    # "omega_slope_surface": (-1.2, 1.2),
    # "omega_surface": (-0.15, 0.15),})[:ndim]

    #rng = ndim * [0.98]
    # rng[0] = (0.0, 0.25)

	corner.corner(ss,
        # s[nskip:, :, :].reshape((-1, ndim)),
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        range=rng, 
        plot_datapoints=True,
        smooth=0.8,
        labels=samples.param_names[:ndim],
    )
	plt.gcf().set_size_inches((20, 20))
	plt.savefig(args.corner)
	plt.clf()

if args.chains is not None:
    fig, axs = plt.subplots(nrows=ndim, ncols=1, figsize=(15, ndim * 6))
    sig1, sig2, sig3 = 0.68, 0.95, 0.9973
    hsig1, hsig2, hsig3 = 0.5 * sig1, 0.5 * sig2, 0.5 * sig3
    for i, p in enumerate(samples.param_names[:ndim]):
        x = np.arange(s[nskip:].shape[0])
        q = np.quantile(
            s[nskip:],
            [
                0.5 - hsig3,
                0.5 - hsig2,
                0.5 - hsig1,
                0.5,
                0.5 + hsig1,
                0.5 + hsig2,
                0.5 + hsig3,
            ],
            axis=1,
        )
        # axs[i].plot(s[:,:,i], lw=0.1, color="tab:blue", label=p, alpha=0.2)
        axs[i].fill_between(x, q[0, :, i], q[6, :, i], color="gray", alpha=0.3, lw=1)
        axs[i].fill_between(x, q[1, :, i], q[5, :, i], color="gray", alpha=0.3, lw=1)
        axs[i].fill_between(x, q[2, :, i], q[4, :, i], color="gray", alpha=0.3, lw=1)
        axs[i].plot(x, q[3, :, i], color="tab:red", lw=2)
        axs[i].set_ylabel(p)
    plt.savefig(args.chains, bbox_inches="tight")
    plt.clf()

if args.logp is not None:
    plt.figure()
    lpp = lp[nskip:, :]
    M = np.max(lpp)
    m = np.median(lpp)
    delta = M - m
    plt.hist(lpp.flat, range=(M - 3 * delta, M), bins=50, density=False)

    # plt.semilogy()
    # lam = np.linspace(M-3*delta, M, 100)
    # lam0 = 0.0
    # plt.plot(lam, 2.0*np.exp(2*(lam-lam0)), lw=2, color="k")

    plt.axvline(m, color="tab:red", lw=2)
    plt.savefig(args.logp, bbox_inches="tight")
    plt.clf()

if args.conv is not None:
    fig, axs = plt.subplots(nrows=ndim, ncols=1, figsize=(15, ndim * 6))
    for i, p in enumerate(samples.param_names[:ndim]):
        ssub = s[nskip:, :, i]
        sig = ssub.std(axis=1).mean()

        # med = np.median(s[:,:,i], axis=1)
        # dev = np.diff(med) * 100.0 / sig
        # axs[i].set_xlabel("sample")

        med = np.median(ssub, axis=1)
        dev = (med - med[-1]) / sig
        axs[i].set_xlabel("last samples")

        axs[i].axhline(0, color="k", ls="--")
        axs[i].plot(dev, color="tab:red")
        axs[i].set_ylim(-0.1, 0.1)
        axs[i].set_ylabel(p)
    plt.savefig(args.conv, bbox_inches="tight")
    plt.clf()

if args.post_csv is not None:
    ibeg = -10
    ssub = s[ibeg:, :, :].reshape((-1, ndim))
    sig1 = 0.68
    hsig1 = 0.5 * sig1
    p_lo, p_med, p_hi = np.quantile(ssub, [0.5 - hsig1, 0.5, 0.5 + hsig1], axis=0)
    p_sig = p_hi - p_lo
    df = pd.DataFrame()
    df["parameter"] = samples.param_names[:ndim]
    df["med"] = p_med
    df["cr68p"] = p_sig
    df = df.set_index("parameter")
    df.to_csv(args.post_csv)
