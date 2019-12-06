#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects solver for single agents
"""

import numpy as np
#from scipy.optimize import fminbound

#from opt_test import build_s_grid, sgrid_on_agrid, get_EVM
from optimizers import v_optimize_single



def v_iter_single(setup,t,EV,female):
    #EVT = np.float32(EV.T)
    
    agrid_s = setup.agrid_s
    sgrid_s = setup.sgrid_s
    
    ind, p = setup.s_ind_s, setup.s_p_s
    
    
    zvals = setup.exogrid.zf_t[t] if female else setup.exogrid.zm_t[t]
    sigma = setup.pars['crra_power']
    beta = setup.pars['beta_t'][t]
    R = setup.pars['R_t'][t]
    
    
    #money = R*agrid_s[:,None] + np.exp(zvals[None,:])
    shp = (agrid_s.size,zvals.size)
    
    
    V_ret, c_opt, s_opt = np.empty_like(shp), np.empty_like(shp), np.empty_like(shp)
    
    
    money_t = (R*agrid_s,np.exp(zvals))
    
    
    V_ret, c_opt, s_opt = v_optimize_single(money_t,sgrid_s,(ind,p,EV),sigma,beta)
    
    def r(x): return x.astype(np.float32)
    
    return r(V_ret), r(c_opt), r(s_opt)#, s_opt/money
