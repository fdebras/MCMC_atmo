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
dates,Vfiles,Ifiles,num_transit,orders,Wmean):

    mod_2D = []
    for mod in reduced_dic["models"]:
        mod_dic = dict(nord=mod.nb,wavelength=mod.Wm,absorption=mod.DD)
        mod_2D_order = prep_func.Model_2D(mod_dic)
        mod_2D_order.create_interp()
        mod_2D.append(mod_2D_order)
        
        
    Rp_Rs = Rp/(Rs)

	### Enter the planet parameters
    P = []
    for i in range(num_transit):
        P_indiv = prep_func.Planet(Rp_Rs,inc,t0,sma,orb_per,w_peri,ecc,limbdark,u_limbdark,dates[i])
	### Generate the transit event
        P_indiv.make_batman()

	### Build weighting window to transform the template of 2D sequence
        P_indiv.make_window()
        P.append(P_indiv)


	# function of interpolation of the model at each phase
    
    # for mod in mod_2D:
    #     mod.build_model(P.window)

    sig_v_int = 1.0   ### Half-width of the window [km/s] (we take half of SPIRou pixel)
    N_v_int   = 7     ### Number of points of the integration
    ddv = np.linspace(-sig_v_int,sig_v_int,N_v_int)
    
    final_model = []
    final_data = []
    for i in range(num_transit):
        tot_indiv = prep_func.total_model(para_dic["Kp"],para_dic["Vsys"],orders[i],Wmean[i],mod_2D,Vfiles[i],Ifiles[i],P[i],ddv)
        tot_indiv.fill_model()
        indiv_model,indiv_data = tot_indiv.bin_model()
        final_model = final_model+indiv_model # concatenate the models to have a list of spectra
        final_data = final_data+indiv_data # same here for the data

    return {
			"data": final_data,
            "model": final_model
        }