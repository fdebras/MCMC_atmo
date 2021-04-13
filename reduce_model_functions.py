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




        
def norm_pts(n_points,Wm,Im,N_bor):
    ### Function to automatically 
    ### 1- Bin the data at n_points
    ### 2- Smooth the resulting data set
    ### 3- Divide the input spectrum by the resulting envelope  
    
    ### Bin the data
    W_bin = []
    I_bin = []
    W_lim = np.linspace(Wm[0],Wm[-1],n_points)
    
    W_bin.append(np.median(Wm[:N_bor]))
    r,cl,cm = stats.sigmaclip(Im[:N_bor],3,3)
    I_bin.append(np.median(r))
    for k in range(n_points-1):
        N_inf = np.argmin(np.abs(Wm - W_lim[k]))
        N_sup = np.argmin(np.abs(Wm - W_lim[k+1]))
        W_bin.append(np.median(Wm[N_inf:N_sup]))
        r,cl,cm = stats.sigmaclip(Im[N_inf:N_sup],3,3)
        I_bin.append(np.median(r))
    W_bin.append(np.median(Wm[-N_bor:]))
    r,cl,cm = stats.sigmaclip(Im[-N_bor:],3,3)
    I_bin.append(np.median(r))
    
    W_bin,I_bin = np.array(W_bin,dtype=float),np.array(I_bin,dtype=float)
    
    ### Interpolate the binned data using a smooth function
    
    #test = interpolate.splrep(W_bin,I_bin,s=0)
    #I_pred_test = interpolate.splev(Wm,test,der=0)
    test = interpolate.interp1d(W_bin,I_bin,kind="linear",fill_value="extrapolate")
    I_pred_test = test(Wm)
    ### Divide by the smooth function
    I_nor  = Im - I_pred_test
    I_nor -= np.max(I_nor)
    return I_nor, I_pred_test, W_bin, I_bin


class reduced_order:
    
    def __init__(self,nb,Wm,Rp,R_s):
        
        self.nb = nb
        self.Wm          = Wm # wavelength model
        self.Rp          = Rp
        self.Vm          = [] # velocity model
        self.DD          = []   ## 1D vector -- Rp as fct of Wm
        
        self.R_s = R_s
        self.models = []
        
        
    def read_adjust_model(self,n_per=5,bin_d=0,n_bin=0):
##### TO IMPROVE


        DF   = -1.0* self.Rp**(2)/self.R_s**(2)

### Get the enveloppe of the distribution
        n_points = n_bin
        N_bor    = 50
		
		### Remove 50% of lowest data pts to get the enveloppe
        I_per = np.percentile(DF,50)
        I_tmp = DF[np.where(DF>I_per)]
        W_tmp = self.Wm[np.where(DF>I_per)]

        
        I_nor, I_pred_test, W_bin, I_bin = norm_pts(n_points,W_tmp,I_tmp,N_bor)
        f = interpolate.interp1d(W_tmp,I_pred_test,kind="linear",fill_value="extrapolate")
        I_fin = f(self.Wm)
        DD  = DF-I_fin#DF*DF/I_fin
        DD -= np.percentile(DD,99)

        self.DD = DD   
        
        
        # c0      = 29979245800.0e-5
        # V_mod   = c0*(self.Wm/self.W_mean-1)
        # self.Vm = V_mod 

        
class reduced:
    
    def __init__(self,model_dic,R_s,orders):

        self.list_ord  = orders  ## List of orders to be used -- after preselection
        self.list_name = []  ## List of the names of the observations
        
        self.Wm = model_dic["freq_nm"]
        self.Rp = model_dic["radius_transm_cm"]
        self.R_s = R_s
        
        
        self.models  = []
        self.N_ord   = 0
        
        self.make_tp = 0

     
            
    def orders_models(self):  
        for i in range(len(self.Wm)):   
            ### First we have to select only the portion of  the model that is of interest for us
            # limits = [self.lambdas[no][0]*0.995,self.lambdas[no][1]*1.005]
            no = self.list_ord[i]
            M  = reduced_order(no,self.Wm[i],self.Rp[i],self.R_s)
            self.models.append(M)
        

        
    def make_all(self,n_per,make_tp=0,bin_d=0,n_bin=0):
        for M in self.models:                        
            M.read_adjust_model(n_per,bin_d,n_bin)
            
    

    




    