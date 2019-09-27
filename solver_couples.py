#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is solver for those who are couples at period 0
"""
import numpy as np

from opt_test import build_s_grid, sgrid_on_agrid, get_EV, v_optimize_multiEV


def vm_period_zero_grid(setup,a0,EV):
    
    agrid = setup.agrid
    sgrid = build_s_grid(agrid,20,0.001,0.1)
    ind, p = sgrid_on_agrid(sgrid,agrid)
    
    
    zf  = setup.exogrid.all_t[0][:,0]
    zm  = setup.exogrid.all_t[0][:,1]
    psi = setup.exogrid.all_t[0][:,2]
    beta = setup.pars['beta']
    sigma = setup.pars['crra_power']

    
    money = a0[:,None] + np.exp(zf[None,:]) + np.exp(zm[None,:])
    
    shp = (setup.na,setup.nexo,setup.ntheta)
    
    
    V_ret, c_opt, s_opt = np.empty(shp), np.empty(shp), np.empty(shp)
    
    umult_vec = setup.u_mult(setup.thetagrid)
    
    for iexo in range(setup.nexo):
        mi = money[:,iexo]
        uadd = psi[iexo]
        
        MEV = np.empty((sgrid.size,setup.ntheta),np.float64)
        
        for itheta in range(setup.ntheta):
            umult = umult_vec[itheta]
            MEV[:,itheta] = (1/umult)*get_EV(ind,p,EV[:,iexo,itheta])
            
        q = v_optimize_multiEV(mi,sgrid,MEV,1,sigma,beta,0)
                
        V_ret[:,iexo,:] = umult_vec[None,:]*q[0] + uadd
        c_opt[:,iexo,:], s_opt[:,iexo,:] = q[1:]
                
    return V_ret, c_opt, s_opt