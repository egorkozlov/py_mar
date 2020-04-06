#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects solver for single agents
"""

import numpy as np
#from scipy.optimize import fminbound

#from opt_test import build_s_grid, sgrid_on_agrid, get_EVM
from optimizers import v_optimize_single, v_optimize_couple, get_EVM



def v_iter_single(setup,t,EV,female,ushift):
    
    agrid_s = setup.agrid_s
    sgrid_s = setup.sgrid_s
    
    
    dtype = setup.dtype
    
    ind, p = setup.vsgrid_s.i, setup.vsgrid_s.wthis
    
    using_pre_u=setup.usinglef_precomputed_u if female else setup.usinglem_precomputed_u
    using_pre_x=setup.usinglef_precomputed_x if female else setup.usinglem_precomputed_x
    zvals = setup.exogrid.zf_t[t] if female else setup.exogrid.zm_t[t]
    ztrend = setup.pars['f_wage_trend'][t] if female else setup.pars['m_wage_trend'][t]
    #sigma = setup.pars['crra_power']
    beta = setup.pars['beta_t'][t]
    R = setup.pars['R_t'][t]
    
    
    #money_t = (R*agrid_s,np.exp(zvals + ztrend))#,np.zeros_like(zvals))    
    #EV_t = (ind,p,EV)
    #V_ret, c_opt, s_opt = v_optimize_single(money_t,sgrid_s,EV_t,sigma,beta,ushift,dtype=dtype)
    
    
    money_t = (R*agrid_s,np.exp(zvals + ztrend),np.zeros_like(zvals))
    ls = np.array([1.0])
    
    
    if EV is None:
        EV = np.zeros((agrid_s.size,zvals.size))
    
    
    V_0, c_opt, x_opt, s_opt, i_opt, _, _ = \
        v_optimize_couple(money_t,sgrid_s,(ind,p,EV[:,:,None,None]),setup.mgrid,
                             using_pre_u[:,None,None],
                             using_pre_x[:,None,None],
                                 ls,beta,ushift,dtype=dtype)
    
    
    
    V_0, c_opt, x_opt, s_opt, i_opt =  \
        (x.squeeze(axis=2) for x in [V_0, c_opt, x_opt, s_opt, i_opt])
    
    
    EVexp = get_EVM(ind,p,EV)
    V_ret = setup.u_single_pub(c_opt,x_opt,ls) + ushift + beta*np.take_along_axis(EVexp,i_opt,0)
    
    def r(x): return x.astype(dtype)
    
    return r(V_ret), r(c_opt), r(x_opt), r(s_opt)
