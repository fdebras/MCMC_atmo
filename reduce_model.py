import numpy as np
from numpy import ma

import sys
import os
import time

from scipy import signal
from scipy import interpolate
from scipy.optimize import curve_fit
import scipy.stats as stats


import astropy
from astropy.io import fits
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting, polynomial


import reduce_model_functions as mod_func
import batman



def reduce(model_dic,R_s,lambdas,orders):


	### Technical parameters for the model
	# Number of percentiles to be set to 0:
	bin_d = 1    ### If bin_d == 1 --> ajust continuum level by binning data and fit with smooth fct
	n_bin = 50
	n_per = 0.1


	lim_data = 1.0
	factor = 0.0005


	MOD = mod_func.reduced(model_dic,R_s,lambdas,orders)
	MOD.list_ord = orders
	MOD.orders_models()

	MOD.make_all(n_per,bin_d,n_bin)

	return {
            "models": MOD.models
        }










        