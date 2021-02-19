import numpy as np
from numpy import ma

import sys
import os
import time

from scipy import interpolate
import batman

import matplotlib.pyplot as plt
import matplotlib

import prepare_likelihood_functions as prep_func



def data_and_model_dict(reduced_dic,para_dic,Rp,Rs,inc,t0,sma,orb_per,w_peri,ecc,limbdark,u_limbdark, \
dates,Vfiles,Ifiles):

    mod_2D = []
    for mod in reduced_dic["models"]:
        mod_dic = dict(velocity=mod.Vm,absorption=mod.DD)
        mod_2D_order = prep_func.Model_2D(mod_dic)
        mod_2D.append(mod_2D_order)
        
    Rp_Rs = Rp/(Rs)

	### Enter the planet parameters
    P = prep_func.Planet(Rp_Rs,inc,t0,sma,orb_per,w_peri,ecc,limbdark,u_limbdark,dates)

	### Generate the transit event
    P.make_batman()

	### Build weighting window to transform the template of 2D sequence
    P.make_window()


	# function of interpolation of the model at each phase
    for mod in mod_2D:
        mod.build_model(P.window)
	
    tot = prep_func.total_model(para_dic["Kp"],para_dic["Vsys"],mod_2D)
    tot.planet = P      ### Planet object previously defined
    
    tot.Vfiles = Vfiles
    tot.Ifiles = Ifiles
    
	# Define integration window (no need to modify this, see note above)
    sig_v_int = 1.0   ### Half-width of the window [km/s] (we take half of SPIRou pixel)
    N_v_int   = 7     ### Number of points of the integration

	### Create instances of the attributes of the CCF object

    tot.ddv    = np.linspace(-sig_v_int,sig_v_int,N_v_int)
    
    tot.models  = mod_2D    ### Model object previously defined

    
    final_model,final_data = tot.bin_model()

    return {
			"data": final_data,
            "model": final_model
        }