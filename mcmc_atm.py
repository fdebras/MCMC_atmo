import os
import sys
import json
import runpy
from argparse import ArgumentParser

import numpy as np
import emcee
import corner

import petit_model

import likelihood_atm as lhood

from mpi4py import MPI
import dill

MPI.pickle.__init__(dill.dumps, dill.loads)

import schwimmbad


def mpi_print(*args):
    if MPI.COMM_WORLD.Get_rank() == 0:
        print(*args)


# Make sure we'll be running libraries on single core
os.environ["OMP_NUM_THREADS"] = "1"


# def make_corr_matrix(args):
#     if args.corr_matrix is None:
#         mpi_print("Using identity correlation matrix")
#         return lhood.make_corr_matrix_identity(args.Jmax, even_only=True)
#     else:
#         mpi_print("Using correlation matrix from", args.corr_matrix)
#         return lhood.make_corr_matrix_from_csv(args.corr_matrix)


likelihood_factories = {
    # "gaussian": lambda args: lhood.FullLikelihood(
    #     lhood.GaussianMassAbundanceLikelihood(),
    #     lhood.MultigaussianMomentsLikelihood(make_corr_matrix(args), args.Jmax),
    # ),
    # "student": lambda args: lhood.FullLikelihood(
    #     lhood.StudentMassAbundanceLikelihood(args.nu),
    #     lhood.MultistudentMomentsLikelihood(make_corr_matrix(args), args.Jmax, args.nu),
    # ),
    "brogi": lambda args: lhood.FullLikelihood(
        lhood.BrogiLikelihood()),
    "Gibson": lambda args: lhood.FullLikelihood(
        lhood.GibsonLikelihood())
}


def parse_cmdline_args():
    parser = ArgumentParser()

    g = parser.add_argument_group("Physical model")
    g.add_argument(
        "--enable-diff-rotation",
        action="store_true",
        help="enable differential rotation",
    )
    
    g = parser.add_argument_group("Statistical model")
    g.add_argument(
        "--data",
        required=True,
        metavar="PY_FILE",
        help="Python file defining planet data (J* observations, ...) to use",
    )
    g.add_argument(
        "--prior",
        required=True,
        metavar="PY_FILE",
        help="Python file defining priors to use",
    )
    g.add_argument(
        "--like",
        choices=likelihood_factories.keys(),
        required=True,
        help="likelihood to use",
    )
    g = parser.add_argument_group("MCMC")
    g.add_argument("--nsteps", "-n", type=int, default=500)
    g.add_argument("--nwalkers", "-w", type=int, default=1000)
    g.add_argument("--out", required=True, help="HDF5 file for MCMC samples output")
    g.add_argument("--seed", type=int, default=42, help="master RNG seed to use")
    g.add_argument(
        "--resume", action="store_true", help="resume from existing MCMC samples output"
    )
    g.add_argument(
        "--abort-exit",
        action="store_true",
        help="use MPI_Abort instead of sys.exit to exit main rank at end",
    )

    g = g.add_mutually_exclusive_group(required=True)
    g.add_argument(
        "--start-range",
        metavar="JSON_FILE",
        help="start walkers from ranges in given JSON file",
    )
    g.add_argument(
        "--start-prior",
        action="store_true",
        help="start walkers from positions drawn from the prior",
    )

    args = parser.parse_args()
    mpi_print(args)
    return args

global args
args = parse_cmdline_args()

# Find out what parameters we're running with ======================================================

param_names = [
        "MMR_H2O",
        # "MMR_CO",
        # "MMR_CO2",
        "Kp",
        "Vsys",
]
        # "kappa_IR",
        # "gamma",
        # "T_int",
        # "T_eq",

if args.enable_diff_rotation:
    param_names += [
        "omega_slope_l0",
        "omega_slope_surface",
        "omega_surface",
    ]
unprior = dict( kappa_IR=0.01,
        gamma=0.012,
        T_int=200,
        T_eq=1000,
        MMR_H2O=0.001,
        MMR_CO=0.0,
        MMR_CO2=0.0,)
ndim = len(param_names)


def unpack_theta(theta):
    assert len(param_names) == len(theta)
    return dict(zip(param_names, theta))


def pack_theta(d):
    return np.array([d[p] for p in param_names])


# Read planet data
planet_data = runpy.run_path(args.data)["make_data"](args)

config = dict(
    p_minbar = -8.0,
    p_maxbar = 2.0,
    n_pressure = 130,
    mass_MJ= planet_data["mass_MJ"],
    radius_RJ=planet_data["radius_RJ"],
    Rs_Rsun = planet_data["Rs_Rsun"],
    gravity_SI= planet_data["gravity_SI"],
    P0_bar = 0.1,
    HHe_ratio=0.275,  # solar ratio
    #now the transit parameters
    inc = planet_data["inc"],
    t0 = planet_data["t0"],
    sma = planet_data["sma"],
    orb_per = planet_data["orb_per"],
    ecc = planet_data["ecc"],
    w_peri = planet_data["w_peri"],
    limbdark = planet_data["limbdark"],
    u_limbdark = planet_data["u_limbdark"],
    dates = planet_data["dates"],
    #the limit of the SPIRou orders
    Wmean = planet_data["Wmean"],
	lambdas = planet_data["lambdas"],
    orders = planet_data["orders"],
    orderstot = planet_data["orderstot"],
	#file with the reduced  datas
	Vfiles = planet_data["Vfiles"],
	Ifiles = planet_data["Ifiles"],
    Stdfiles = planet_data["Stdfiles"],
    num_transit = planet_data["num_transit"],
)
config.update(unprior)
# Posterior and model ==============================================================================


class Posterior:
    def __init__(
        self, prior, like, model, planet_data, config, enable_diff_rotation=False
    ):
        self.prior = prior
        self.like = like
        self.model = model
        self.planet_data = planet_data
        self.config = config
        self.enable_diff_rotation = enable_diff_rotation

    def ln_post_vector(self, theta):
        return self.ln_post(unpack_theta(theta))

    def ln_post(self, d):
        ln_p = self.prior.ln_prior(d)
        if not np.isfinite(ln_p):
            return -np.inf
        else:
            return ln_p + self.ln_like(d)

    def ln_like(self, d):
        d = dict(d)
        try:
            resid = self.model.prepare_likelihood(d)
        except RuntimeError:
            return -np.inf
        self.like.corr = resid
        return self.like.ln_like()


def make_posterior(config):
    model = petit_model.Model(config)

    # Import prior module and execute it
    prior = runpy.run_path(args.prior)["make_prior"](args)

    assert tuple(sorted(param_names)) == tuple(
        sorted(prior.params)
    ), "params mismatch between model: {} and priors: {}".format(
        str(param_names), str(prior.params)
    )

    likelihood = likelihood_factories[args.like](args)
    return Posterior(
        prior,
        likelihood,
        model,
        planet_data,
        config,
        enable_diff_rotation=args.enable_diff_rotation,
    )


# NOTE: `post` (more specifically, `model`) MUST be _constructed_ on each MPI rank,
# otherwise MPIPool() will pickle and share instances of it across ranks which breaks.
post = make_posterior(config)

# main() will only be called by MPI main rank, which will dispatch work to worker ranks
# So everything in main() is serial, unless explicitly dispatched to the pool
def main(pool):

    # Walker initialization ============================================================================

    if args.start_range is not None:
        with open(args.start_range, "r") as f:
            param_ranges = json.load(f)

        random_params0 = lambda: {
            p_name: np.random.uniform(low=p_low, high=p_high)
            for p_name, (p_low, p_fid, p_high) in param_ranges.items()
        }
    else:
        assert args.start_prior
        random_params0 = lambda: post.prior.rvs()

    def make_pos0(i):
        np.random.seed(args.seed + 1664 + i)
        max_tries = 1000
        for _ in range(max_tries):
            params = random_params0()
            # print(params)
            if np.isfinite(post.ln_post(params)):
                print("Initial position for walker {:05d} = {}".format(i, str(params)))
                return pack_theta(params)
        raise ValueError(
            "cannot find a start position with finite posterior"
            f" for walker #{i} after {max_tries} tries, check your start position options"
        )

    assert args.out.endswith(".h5")

    if not args.resume:
        # Need to generate initial seed positions
        print("Will generate walker starting positions...")
        pos0 = pool.map(make_pos0, range(args.nwalkers))
        np.save(args.out.replace(".h5", "") + "_pos0.npy", pos0)
    else:
        # Resuming from previous run
        print("Will attempt resuming from previous MCMC samples in", args.out)
        # Set initial positions to None for sampler to reuse from backend
        pos0 = None

    print("Initial walker positions generated")

    # Backend ==========================================================================================

    np.random.seed(args.seed)

    # Set up the backend
    # Don't forget to clear it in case the file already exists
    backend = emcee.backends.HDFBackend(args.out)

    if not args.resume:
        # Make sure to reset the backend to overwrite any existing output if we're not resuming
        backend.reset(args.nwalkers, ndim)

    # Start sampling ===================================================================================

    def ln_post(theta):
        return post.ln_post_vector(theta)

    sampler = emcee.EnsembleSampler(
        args.nwalkers, ndim, ln_post, pool=pool, backend=backend
    )
    print("Starting sampling...")
    sampler.run_mcmc(pos0, args.nsteps, progress=True)
    print("Done sampling!")


if __name__ == "__main__":
    pool = schwimmbad.MPIPool()
    if pool.is_master():
        # Start processing on master rank only
        main(pool)
        print("Exited main() on main rank")
        if args.abort_exit:
            print("abort_exit requested, will now exit with MPI_Abort()")
            os.sync()
            MPI.COMM_WORLD.Abort()
        else:
            sys.exit(0)
    else:
        # Other ranks wait for MPIPool() to dispatch work to them
        pool.wait()
        sys.exit(0)
