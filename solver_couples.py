#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is solver for those who are couples at period 0
"""
import numpy as np
from numba import jit

from opt_test import build_s_grid, sgrid_on_agrid, get_EV, get_EVM

#from opt_test import v_optimize_multiEV as v_optimize
from platform import system

if system() != 'Darwin':
    from opt_test import v_optimize_MEV_cp as v_optimize
    from opt_test import v_optimize_MEV_cp_massive as v_optimizeM
else:
    from opt_test import v_optimize_MEV_np as v_optimize
    from opt_test import v_optimize_MEV_np_massive as v_optimizeM
    
    
#from aot_test import vopt_MEV

#@jit(nopython=True)
def vm_period_zero_grid(setup,a0,EV):
    
    agrid = setup.agrid
    sgrid = build_s_grid(agrid,10,0.001,0.1)
    ind, p = sgrid_on_agrid(sgrid,agrid)
    
    
    zf  = setup.exogrid.all_t[0][:,0]
    zm  = setup.exogrid.all_t[0][:,1]
    psi = setup.exogrid.all_t[0][:,2]
    beta = setup.pars['beta']
    sigma = setup.pars['crra_power']

    
    money = a0[:,None] + np.exp(zf[None,:]) + np.exp(zm[None,:])
    
    shp = (setup.na,setup.nexo,setup.ntheta)
    
    money,sgrid,EV,sigma,beta = (np.float32(x) for x in (money,sgrid,EV,sigma,beta))
    
    V_ret, c_opt, s_opt = np.empty(shp), np.empty(shp), np.empty(shp)
    
    umult_vec = setup.u_mult(setup.thetagrid)
    
    MMEV = (1/umult_vec[None,None,:])*get_EVM(ind,p,EV)
    
    
    
    
    nbatch = 20
    istart = 0
    ifinish = nbatch
    
    
    
    for ibatch in range(int(np.floor(setup.nexo/nbatch))):
        money_i = money[:,istart:ifinish]
        MMEV_i = MMEV[:,istart:ifinish]
        V_pure_i, c_opt_i, s_opt_i = v_optimizeM(money_i,sgrid,MMEV_i,sigma,beta)
        V_ret_i = umult_vec[None,None,:]*V_pure_i + psi[None,istart:ifinish,None]
        
        
        V_ret[:,istart:ifinish,:] = V_ret_i
        c_opt[:,istart:ifinish,:] = c_opt_i
        s_opt[:,istart:ifinish,:] = s_opt_i
        
        istart = ifinish
        ifinish = ifinish+nbatch if ifinish+nbatch < setup.nexo else setup.nexo
        
    
    
    
    #V_pure, c_opt, s_opt = v_optimizeM(money,sgrid,MMEV,sigma,beta)
    
    #V_ret = umult_vec[None,None,:]*V_pure + psi[None,:,None]
    
           
    return V_ret, c_opt, s_opt



def vm_period_zero_grid_loop(setup,a0,EV):
    
    agrid = setup.agrid
    sgrid = build_s_grid(agrid,10,0.001,0.1)
    ind, p = sgrid_on_agrid(sgrid,agrid)
    
    
    zf  = setup.exogrid.all_t[0][:,0]
    zm  = setup.exogrid.all_t[0][:,1]
    psi = setup.exogrid.all_t[0][:,2]
    beta = setup.pars['beta']
    sigma = setup.pars['crra_power']

    
    money = a0[:,None] + np.exp(zf[None,:]) + np.exp(zm[None,:])
    
    shp = (setup.na,setup.nexo,setup.ntheta)
    
    money,sgrid,EV,sigma,beta = (np.float32(x) for x in (money,sgrid,EV,sigma,beta))
    
    V_ret, c_opt, s_opt = np.empty(shp), np.empty(shp), np.empty(shp)
    
    umult_vec = setup.u_mult(setup.thetagrid)
    
    for iexo in range(setup.nexo):
        mi = money[:,iexo]
        uadd = psi[iexo]       
        
        MEV = (1/umult_vec[None,:])*get_EVM(ind,p,EV[:,iexo,:])        
        
        q = v_optimize(mi,sgrid,MEV,np.float32(1.0),sigma,beta,np.float32(0.0))
        
                
        V_ret[:,iexo,:] = umult_vec[None,:]*q[0] + uadd
        c_opt[:,iexo,:], s_opt[:,iexo,:] = q[1:]
    
           
    return V_ret, c_opt, s_opt