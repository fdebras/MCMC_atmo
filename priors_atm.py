""" Main module for MCMC priors

All priors are objects which should expose:

    - A `params` attribute, which holds the list of parameter names that this prior handles
    - A `ln_prior()` method, which evaluates the log-prior on a given dictionary of parameter values
    - A `rvs()` method, which draws a random variate sample according to the prior (returns a dict)
"""

import numpy as np
from scipy import stats


class FullPrior(object):
    def __init__(self, priors):
        self.priors = priors
        all_params = tuple(sum([prior.params for prior in priors], []))
        assert len(all_params) == len(
            set(all_params)
        ), "sub-priors of FullPrior should not have overlapping parameters"
        self.params = all_params


    def ln_prior(self, param_dict):
        ln_p = 0.0
        for prior in self.priors:
            ln_p_this_prior = prior.ln_prior(param_dict)
            if np.isfinite(ln_p_this_prior):
                ln_p += ln_p_this_prior
            else:
                return -np.inf
        return ln_p

    def rvs(self):
        d = {}
        for prior in self.priors:
            d.update(prior.rvs())
        return d


def ln_prior_beta(x, a, b):
    if not (0 < x < 1):
        return -np.inf
    else:
        return (a - 1.0) * np.log(x) + (b - 1.0) * np.log(1.0 - x)


class SimplexZonesPrior(object):
    params = [
        "core_fraction",
        "lambda_diffuse",
        "lambda_min_layered",
        "lambda_max_layered",
    ]

    def __init__(self, a=18.0, b=8.0):
        self.a, self.b = a, b
        self.beta = stats.beta(a=a, b=b)

    def _is_inside_simplex(self, d):
        return (
            0.0
            < d["core_fraction"]
            < d["lambda_diffuse"]
            < d["lambda_min_layered"]
            < d["lambda_max_layered"]
            < 1.0
        )

    def ln_prior(self, param_dict):
        if not self._is_inside_simplex(param_dict):
            return -np.inf
        else:
            ln_p_lambda_max = ln_prior_beta(
                param_dict["lambda_max_layered"], a=self.a, b=self.b
            )
            # Careful about volume of simplex: depends on lambda_max_layered
            # Because each variable is uniform between 0 and lambda_max_layered
            ln_p_simplex_cond_lambda_max = -3.0 * np.log(
                param_dict["lambda_max_layered"]
            )
            return ln_p_lambda_max + ln_p_simplex_cond_lambda_max

    def rvs(self):
        # First draw lambda_max_layered
        lam_max = self.beta.rvs()
        # Then get the other ones with rejection sampling
        while True:
            d = {
                "core_fraction": np.random.uniform(low=0, high=lam_max),
                "lambda_diffuse": np.random.uniform(low=0, high=lam_max),
                "lambda_min_layered": np.random.uniform(low=0, high=lam_max),
                "lambda_max_layered": lam_max,
            }
            if self._is_inside_simplex(d):
                return d


class SimplexUniformT_Priors(object):
    params = ["T_eq",
#        "kappa_IR",
#        "gamma",
#        "T_int",
#        "T_eq",
    ]

    def __init__(self, T_eq_min,T_eq_max):#kappa_min=-5.0, kappa_max=-2.0, gamma_min=-2.0,gamma_max=2.0,T_int_min=0.0,T_int_max=500.0, \
#    T_eq_min=800.0,T_eq_max=2500.0):
    
#        self.kappa_min = kappa_min
#        self.kappa_max = kappa_max
        
#        self.gamma_min = gamma_min
#        self.gamma_max = gamma_max

        
#        self.T_int_min = T_int_min
#        self.T_int_max = T_int_max
        
        self.T_eq_min = T_eq_min
        self.T_eq_max = T_eq_max

    

    def _is_inside_simplex(self, d):
        return (  self.T_eq_min<d["T_eq"]<self.T_eq_max
#self.kappa_min<d["kappa_IR"]<self.kappa_max
 #               and self.gamma_min<d["gamma"]<self.gamma_max
  #              and self.T_int_min<d["T_int"]<self.T_int_max
   #             and self.T_eq_min<d["T_eq"]<self.T_eq_max
        )

    def ln_prior(self, param_dict):
        if self._is_inside_simplex(param_dict):
            return 0.0
        else:
            return -np.inf

    def rvs(self):
        while True:
#            kir = np.random.uniform(low=self.kappa_min, high=self.kappa_max)
            
#            gamma = np.random.uniform(low=self.gamma_min, high=self.gamma_max)
            d = {"T_eq": np.random.uniform(low=self.T_eq_min, high=self.T_eq_max),
#                "kappa_IR": kir,
#                "gamma": gamma,
#                "T_int": np.random.uniform(low=self.T_int_min, high=self.T_int_max),
#                "T_eq": np.random.uniform(low=self.T_eq_min, high=self.T_eq_max),
            }
            if self._is_inside_simplex(d):
                return d


class SimplexUniformAbund_Priors(object):
    params = [
        "MMR_H2O",
        "MMR_CO",
        # "MMR_CO2",
    ]

    def __init__(self, H2O_min, H2O_max,CO_min,CO_max):
    # , CO2_min=-8.0, CO2_max=-2.0,CO_min=-8.0, CO_max=-2.0):

        self.H2O_min = H2O_min
        self.H2O_max = H2O_max
        
        # self.CO2_min = CO2_min
        # self.CO2_max = CO2_max
        
        self.CO_min = CO_min
        self.CO_max = CO_max


    def _is_inside_simplex(self, d):

        return (self.H2O_min<d["MMR_H2O"]<self.H2O_max
                # and self.CO2_min<d["MMR_CO2"]<self.CO2_max
                and self.CO_min<d["MMR_CO"]<self.CO_max
        )

    def ln_prior(self, param_dict):
        if self._is_inside_simplex(param_dict):
            return 0.0
        else:
            return -np.inf

    def rvs(self):
        while True:
            H2O=  np.random.uniform(low=self.H2O_min, high=self.H2O_max)
            
            CO= np.random.uniform(low=self.CO_min, high=self.CO_max)
            
            # CO2= np.random.uniform(low=self.CO2_min, high=self.CO2_max)        
            d = {
                "MMR_H2O": H2O,
                "MMR_CO": CO,
                # "MMR_CO2": CO2,
            }
            if self._is_inside_simplex(d):
                return d

class SimplexUniformTransit_Priors(object):
    params = [
        "Kp",
        "Vsys"
    ]

    def __init__(self, Kp_min,Kp_max, \
    Vsys_min,Vsys_max):
        self.Kp_min = Kp_min
        self.Kp_max = Kp_max
        
        self.Vsys_min = Vsys_min
        self.Vsys_max = Vsys_max

    

    def _is_inside_simplex(self, d):

        return (self.Kp_min<d["Kp"]<self.Kp_max
                and self.Vsys_min<d["Vsys"]<self.Vsys_max
        )

    def ln_prior(self, param_dict):
        if self._is_inside_simplex(param_dict):
            return 0.0
        else:
            return -np.inf

    def rvs(self):
        while True:
            d = {
                "Kp": np.random.uniform(low=self.Kp_min, high=self.Kp_max),
                "Vsys": np.random.uniform(low=self.Vsys_min, high=self.Vsys_max),
            }
            if self._is_inside_simplex(d):
                return d



class EntropyJumpExpPrior(object):
    params = ["delta_S_layered"]

    def __init__(self, delta_S_exp_scale=3.0e-4, delta_S_cutoff=3e-1):
        self.delta_S_exp_scale = delta_S_exp_scale
        self.delta_S_cutoff = delta_S_cutoff
        self.expon = stats.expon(scale=self.delta_S_exp_scale)

    def ln_prior(self, param_dict):
        if not (0.0 < param_dict["delta_S_layered"] < self.delta_S_cutoff):
            return -np.inf
        else:
            return -param_dict["delta_S_layered"] / self.delta_S_exp_scale

    def rvs(self):
        while True:
            delta_S = self.expon.rvs()
            if delta_S < self.delta_S_cutoff:
                return {"delta_S_layered": delta_S}


class DifferentialRotationPrior(object):
    params = ["omega_slope_l0", "omega_slope_surface", "omega_surface"]

    def __init__(self, omega_surface_absmax=0.2):
        assert omega_surface_absmax > 0.0
        self.omega_surface_absmax = omega_surface_absmax

    def ln_prior(self, param_dict):
        if not (np.abs(param_dict["omega_surface"]) < self.omega_surface_absmax):
            return -np.inf
        else:
            return 0.0

    def rvs(self):
        return {
            "omega_surface": np.random.uniform(
                low=-self.omega_surface_absmax, high=self.omega_surface_absmax
            ),
            "omega_slope_l0": np.random.normal(),
            "omega_slope_surface": np.random.normal(),
        }
