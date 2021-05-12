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

import batman
import matplotlib.pyplot as plt





class reduced_order:
    
    def __init__(self,nb,Wm,Rp,R_s):
        
        self.nb = nb
        self.Wm          = Wm # wavelength model
        self.Rp          = Rp
        self.DD          = []   ## 1D vector -- Rp as fct of Wm
        
        self.R_s = R_s
        self.models = []
        
        
    def read_adjust_model(self,n_per,bin_d,n_bin):
##### TO IMPROVE


        DF   = -1.0* self.Rp**(2)/self.R_s**(2)


        DD = DF
        DD -= np.percentile(DD,99)


        self.DD = DD   
        
        
        # c0      = 29979245800.0e-5
        # V_mod   = c0*(self.Wm/self.W_mean-1)
        # self.Vm = V_mod 

        
class reduced:
    
    def __init__(self,model_dic,R_s,orderstot):

        self.list_ord  = orderstot  ## List of orders to be used -- after preselection
        self.list_name = []  ## List of the names of the observations
        
        self.Wm = model_dic["freq_nm"]
        self.Rp = model_dic["radius_transm_cm"]
        self.R_s = R_s
        
        
        self.models  = []
        self.N_ord   = 0


     
            
    def orders_models(self):  
        for i in range(len(self.list_ord)):   
            ### First we have to select only the portion of  the model that is of interest for us
            # limits = [self.lambdas[no][0]*0.995,self.lambdas[no][1]*1.005]
            no = self.list_ord[i]
            M  = reduced_order(no,self.Wm[i],self.Rp[i],self.R_s)
            self.models.append(M)
        

        
    def make_all(self,n_per,bin_d,n_bin):
        for M in self.models:                        
            M.read_adjust_model(n_per,bin_d,n_bin)
            
    

    




    
