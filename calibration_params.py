#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  8 13:11:16 2020

@author: egorkozlov
"""

import numpy as np
from collections import OrderedDict

def calibration_params(xin=None,xfix=None):
    # NOTE: xfix overwrites xin
    
    
    # format is name: (lb,ub,xinit)
    # I am not sure if this should be ordered or not but let's do ordered
    # just in case...
    
    params = OrderedDict(
              akept=(0.65, 0.9, 0.1),
              sigma_psi=(0.45, 1.0, 0.01),#sigma_psi=(0.005, 0.8, 0.01),
              sigma_psi_mult=(1.60, 3.5, 0.2),#sigma_psi_mult=(0.5, 5.0, 0.02),
              pmeet=(0.33, 0.64, 0.01),#pmeet=(0.1, 1.0, 0.7),
              util_alp=(2.8504, 6.36, 0.01),#util_alp=(0.01, 0.4, 0.25),
              u_shift_mar = (0.0212377, 0.212379, 0.001),
              z_drift = (-0.129, -0.04, -0.001),#z_drift = (-0.3, 0.0, -0.1)
              util_xi=(1.05, 1.22, 0.001)
                        )
                        
    #params = OrderedDict(
     #         akept=(0.907274, 0.907276, 0.1),
      #        sigma_psi=(0.671454, 0.671456, 0.01),#sigma_psi=(0.005, 0.8, 0.01),
       #       sigma_psi_mult=(1.6691, 1.6693, 0.2),#sigma_psi_mult=(0.5, 5.0, 0.02),
        #      pmeet=(0.25998, 0.26, 0.01),#pmeet=(0.1, 1.0, 0.7),
         #     util_alp=(0.606655, 0.606657, 0.01),#util_alp=(0.01, 0.4, 0.25),
          #    u_shift_mar = (0.0212377, 0.0212379, 0.001),#u_shift_mar = (0.0, 0.5, 0.0001),
           #   z_drift = (-0.0506905,-0.0506903, -0.001),#z_drift = (-0.3, 0.0, -0.1)
            #  util_xi=(1.07131,1.07133, 0.001)
             #           )
                
                
                        
   
    # update params is some dict with values was supplied
    if xin is not None:
        assert type(xin) is dict
        for key in xin:
            if type(xin[key]) is tuple:
                params[key] = xin[key]
            else:
                assert key in params
                v_old = params[key]
                params[key] = (v_old[0],v_old[1],xin[key])
                
    if xfix is not None:
        assert type(xfix) is dict
        for key in xfix:
            assert key in params
            if xin is not None and key in xin and xin[key] != xfix[key]:
                print('Warning: xfix overwrites xin')
            xval = xfix[key]
            params[key] = (xval,xval,xval)
                           
    
    keys_fixed, x_fixed = list(), list()
    
    keys, lb, ub, x0 = list(), list(), list(), list()
    
    for x in params:    
        lb_here = params[x][0]
        ub_here = params[x][1]
        x_here = params[x][2]
        
        if np.abs(ub_here - lb_here) < 1e-8: # if ub and lb are equal
            assert lb_here <= x_here <= ub_here
            keys_fixed.append(x)
            x_fixed.append(x_here)
        else:        
            keys.append(x)
            lb.append(params[x][0])
            ub.append(params[x][1])
            x0.append(params[x][2])
                      
            
    lb, ub, x_here, x_fixed = (np.array(x) for x in (lb,ub,x_here,x_fixed))
    
    def translator(x):
        # in case x has a weird dimensionality
        try:
            x.squeeze()
        except:
            pass
        
        assert len(keys) == len(x), 'Wrong lenght of x!'
        d_var = dict(zip(keys,x))
        if len(keys_fixed) > 0:
            d_fixed = dict(zip(keys_fixed,x_fixed))
            d_var.update(d_fixed)
        return d_var
    
    return lb, ub, x0, keys, translator
    