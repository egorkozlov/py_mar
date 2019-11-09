#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects solver for single agents
"""

import numpy as np
#from scipy.optimize import fminbound

#from opt_test import build_s_grid, sgrid_on_agrid, get_EVM
from optimizers import v_optimize



def v_iter_single(setup,EV,female):
    #EVT = np.float32(EV.T)
    
    agrid = setup.agrid
    sgrid = setup.sgrid
    
    ind, p = setup.s_ind, setup.s_p
    
    
    zvals = setup.exogrid.zf_t[0] if female else setup.exogrid.zm_t[0]
    sigma = setup.pars['crra_power']
    beta = setup.pars['beta']
    R = setup.pars['R']
    
    
    money = R*agrid[:,None] + np.exp(zvals[None,:])
    shp = (agrid.size,zvals.size)
    
    
    V_ret, c_opt, s_opt = np.empty_like(shp), np.empty_like(shp), np.empty_like(shp)
    
    
    money_t = (R*agrid,np.exp(zvals))
    
    
    V_ret, c_opt, s_opt = v_optimize(money_t,sgrid,(ind,p,EV),sigma,beta)
    
    
    
    return V_ret, s_opt, s_opt/money
