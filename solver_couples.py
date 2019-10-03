#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is solver for those who are couples at period 0
"""
import numpy as np
from timeit import default_timer


from optimizers import build_s_grid, sgrid_on_agrid, get_EVM
from optimizers import v_optimize as v_optimizeM

from platform import system

if system() != 'Darwin':    
    nbatch_def = 200
else:
    nbatch_def = 1


#@jit(nopython=True)
def vm_period_zero_grid_massive(setup,a0,EV,nbatch=nbatch_def,verbose=False):
    
    if verbose: start = default_timer()
    
    agrid = setup.agrid
    sgrid = build_s_grid(agrid,10,0.001,0.1)
    ind, p = sgrid_on_agrid(sgrid,agrid)
    
    
    zf  = setup.exogrid.all_t[0][:,0]
    zm  = setup.exogrid.all_t[0][:,1]
    psi = setup.exogrid.all_t[0][:,2]
    beta = setup.pars['beta']
    sigma = setup.pars['crra_power']

    labor_income = np.exp(zf) + np.exp(zm)
    
    money = a0[:,None] + labor_income[None,:]
    
    shp = (setup.na,setup.nexo,setup.ntheta)
    
    # type conversion to keep everything float32
    money,sgrid,EV,sigma,beta = (np.float32(x) for x in (money,sgrid,EV,sigma,beta))
    
    V_ret, c_opt, s_opt = np.empty(shp), np.empty(shp), np.empty(shp)
    
    umult_vec = setup.u_mult(setup.thetagrid)
    
    #MMEV = (1/umult_vec[None,None,:])*get_EVM(ind,p,EV)
    # this is required 
    
    EV_resc = (1/umult_vec[None,None,:])*EV
    
    if verbose: print('MMEV computed after {} sec'.format(default_timer()-start))
    
    istart = 0
    ifinish = nbatch if nbatch < setup.nexo else setup.nexo
    
    # this natually splits everything onto slices
    
    for ibatch in range(int(np.floor(setup.nexo/nbatch))):
        #money_i = money[:,istart:ifinish]
        
        money_t = (a0,labor_income[istart:ifinish])
        EV_t = (ind,p,EV_resc[:,istart:ifinish,:])
        
        #MMEV_i = MMEV[:,istart:ifinish]
        #V_pure_i, c_opt_i, s_opt_i = v_optimizeM(money_i,sgrid,MMEV_i,sigma,beta)
        V_pure_i, c_opt_i, s_opt_i = v_optimizeM(money_t,sgrid,EV_t,sigma,beta)
        V_ret_i = umult_vec[None,None,:]*V_pure_i + psi[None,istart:ifinish,None]
        
        
        V_ret[:,istart:ifinish,:] = V_ret_i
        c_opt[:,istart:ifinish,:] = c_opt_i
        s_opt[:,istart:ifinish,:] = s_opt_i
        
        istart = ifinish
        ifinish = ifinish+nbatch if ifinish+nbatch < setup.nexo else setup.nexo
        
        if verbose: print('Batch {} done at {} sec'.format(ibatch,default_timer()-start))
    
    return V_ret, c_opt, s_opt


# this is equivalent to the above function with nbatch = 1
'''
def vm_period_zero_grid_loop(setup,a0,EV):
    
    
    if system() != 'Darwin':
        from opt_test import v_optimize_MEV_cp as v_optimize
    else:
        from opt_test import v_optimize_MEV_np as v_optimize
    
    
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
        
        q = v_optimizeM(mi,sgrid,MEV,sigma,beta)
        #q = v_optimize(mi,sgrid,MEV,sigma,beta)
        
        
                
        V_ret[:,iexo,:] = umult_vec[None,:]*q[0] + uadd
        c_opt[:,iexo,:], s_opt[:,iexo,:] = q[1:]
    
           
    return V_ret, c_opt, s_opt
'''