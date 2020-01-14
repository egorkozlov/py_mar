#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is solver for those who are couples at period 0
"""
import numpy as np
from timeit import default_timer


from optimizers import get_EVM
from optimizers import v_optimize_couple

from platform import system

if system() != 'Darwin' and system() != 'Windows':    
    nbatch_def = 500
    use_cp = True
    
elif system() == 'Windows':
    
    nbatch_def = 17
    use_cp = True
    
else:
    
    nbatch_def = 17
    use_cp = False

def v_iter_couple(setup,t,EV_tuple,ushift,nbatch=nbatch_def,verbose=False):
    
    if verbose: start = default_timer()
    
    agrid = setup.agrid_c
    sgrid = setup.sgrid_c
    ind, p = setup.vsgrid_c.i, setup.vsgrid_c.wthis
    
    
    EV_by_l, EV_fem_by_l, EV_mal_by_l = EV_tuple    
    
    ls = setup.ls_levels
    us = setup.ls_utilities
    #pd = setup.ls_pdown
    nls = len(ls)
    
    
    # type conversion is here
    
    zf  = setup.exogrid.all_t[t][:,0]
    zm  = setup.exogrid.all_t[t][:,1]
    zftrend = setup.pars['f_wage_trend'][t]
    zmtrend = setup.pars['m_wage_trend'][t]

    psi = setup.exogrid.all_t[t][:,2]
    beta = setup.pars['beta_t'][t]
    sigma = setup.pars['crra_power']
    R = setup.pars['R_t'][t]

    
    
    wf = np.exp(zf + zftrend)
    wm = np.exp(zm + zmtrend)
    
    #labor_income = np.exp(zf) + np.exp(zm)
    
    #money = R*agrid[:,None] + wf[None,:] 
    
    
    nexo = setup.pars['nexo_t'][t]
    shp = (setup.na,nexo,setup.ntheta)
    
    # type conversion to keep everything float32
    sgrid,sigma,beta = (np.float32(x) for x in (sgrid,sigma,beta))
    
    V_couple, c_opt, s_opt = np.empty((3,)+shp,np.float32)
    i_opt, il_opt = np.empty(shp,np.int16), np.empty(shp,np.int16)
    
    V_all_l = np.empty(shp+(nls,),dtype=np.float32)
    
    theta_val = np.float32(setup.thetagrid)
    umult_vec = setup.u_mult(theta_val)
    
    # the original problem is max{umult*u(c) + beta*EV}
    # we need to rescale the problem to max{u(c) + beta*EV_resc}
    
    istart = 0
    ifinish = nbatch if nbatch < nexo else nexo
    
    # this natually splits everything onto slices
    
    for ibatch in range(int(np.ceil(nexo/nbatch))):
        #money_i = money[:,istart:ifinish]
        assert ifinish > istart
        
        money_t = (R*agrid, wf[istart:ifinish], wm[istart:ifinish])
        EV_t = (ind,p,EV_by_l[:,istart:ifinish,:,:])
        
        V_pure_i, c_opt_i, s_opt_i, i_opt_i, il_opt_i, V_all_l_i = \
           v_optimize_couple(money_t,sgrid,umult_vec,EV_t,sigma,beta,ls,us,ushift)
        V_ret_i = V_pure_i + psi[None,istart:ifinish,None]
        
        
        
        
        
        V_couple[:,istart:ifinish,:] = V_ret_i
        c_opt[:,istart:ifinish,:] = c_opt_i
        s_opt[:,istart:ifinish,:] = s_opt_i
        i_opt[:,istart:ifinish,:] = i_opt_i
        il_opt[:,istart:ifinish,:] = il_opt_i
        V_all_l[:,istart:ifinish,:,:] = V_all_l_i
        
        
        istart = ifinish
        ifinish = ifinish+nbatch if ifinish+nbatch < nexo else nexo
        
        if verbose: print('Batch {} done at {} sec'.format(ibatch,default_timer()-start))
    
    
    assert np.all(c_opt > 0)
    
    psi_r = psi[None,:,None].astype(np.float32)
    
    # finally obtain value functions of partners
    uf, um = setup.u_part(c_opt,theta_val[None,None,:])
    uc = setup.u_couple(c_opt,theta_val[None,None,:])
    u_pub = us[il_opt]
    
    
    EVf_all, EVm_all, EV_all  = (get_EVM(ind,p,x) for x in (EV_fem_by_l, EV_mal_by_l,EV_by_l))
    V_fem = uf + psi_r + u_pub + ushift+ beta*np.take_along_axis(np.take_along_axis(EVf_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    V_mal = um + psi_r + u_pub + ushift+ beta*np.take_along_axis(np.take_along_axis(EVm_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    V_all = uc + psi_r + u_pub + ushift+ beta*np.take_along_axis(np.take_along_axis(EV_all,i_opt[...,None],0),il_opt[...,None],3).squeeze(axis=3)
    def r(x): return x.astype(np.float32)
    
    assert np.allclose(V_all,V_couple,atol=1e-5,rtol=1e-4)
    
    return r(V_couple), r(V_fem), r(V_mal), r(c_opt), r(s_opt), il_opt, r(V_all_l)




