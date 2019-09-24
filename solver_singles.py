#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects solver for single agents
"""

import numpy as np
from scipy.optimize import fminbound


def v_period_zero_grid(setup,a0,EV,female):
    # this takes gender as argument so should be called twice
    
    agrid = setup.agrid
    zvals = setup.exogrid.zf_t[0] if female else setup.exogrid.zm_t[0]
    beta = setup.pars['beta']
    def u(c): return setup.u(c)
    
    
    income = a0[:,None] + np.exp(zvals[None,:])
    
    
    def neg_total_u(s,inc,EVg):
        c = inc - s
        assert np.all(c > 0)
        EV = np.interp(s,agrid,EVg)        
        return -(u(c) + beta*EV)
    
    
    smin = np.zeros_like(income)
    smax = 0.9*income
    
    EVT = EV.T
    
    s_opt = np.array([
                         [
                          fminbound( lambda x : neg_total_u(x,income[j,i],EVval),smin[j,i],smax[j,i] ) 
                          for i, EVval in enumerate(EVT)
                         ]
                       for j in range(a0.size)
                      ])
    
    
    V_ret =  np.array([-neg_total_u(s_opt[:,i],income[:,i],EV) for i, EV in enumerate(EVT)]).T
    
    return V_ret, s_opt, s_opt/income