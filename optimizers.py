#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects selected routines from opt_test.py

The aim is keep each function readable, and either njit-able or using cuda

Both versions are needed as it is not clear whether we'll actually use cuda

"""

import numpy as np
from numba import jit#, prange, cuda, float32
from platform import system
from aux_routines import cp_take_along_axis

#if system() != 'Darwin' and system() != 'Windows':
if system() != 'Darwin':
    import cupy as cp
    ucp = True
else:
    ucp = False


#@jit('float64[:](float64[:], int64, float64, float64)',nopython=True)
def build_s_grid(agrid,n_between,da_min,da_max):
    sgrid = np.array([0.0],np.float64)
    for j in range(agrid.size-1):
        step = (agrid[j+1] - agrid[j])/n_between
        if step >= da_min and step <= da_max:
            s_add = np.linspace(agrid[j],agrid[j+1],n_between)[:-1]
        elif step < da_min:
            s_add = np.arange(agrid[j],agrid[j+1],da_min)
        elif step > da_max:
            s_add = np.arange(agrid[j],agrid[j+1],da_max)
        sgrid = np.concatenate((sgrid,s_add))
    
    sgrid = np.concatenate((sgrid,np.array([agrid[-1]])))
            
    if sgrid[0] == sgrid[1]: 
        sgrid = sgrid[1:]
        
    return sgrid


def get_EVM(ind,wthis,EVin,use_cp=False):
    # this essenitally doubles VecOnGrid.apply method
    # but we cannot deprecate it unless we are sure VecOnGrid works on GPU 
    # correctly
    mr = cp if use_cp else np
    
    ev_aux_shape  = EVin.shape[1:]
    shap = (ind.size,) + ev_aux_shape
    EVout = mr.empty(shap,mr.float32)
    
    pb = wthis.reshape(((wthis.size,)+(1,)*len(ev_aux_shape)))
    EVout[:] = pb*EVin[ind,...] + (1-pb)*EVin[ind+1,...]    
    
    return EVout


def v_optimize_couple(money,sgrid,umult,EV,sigma,beta,ls,us,ushift,use_cp=ucp,compare=False,force_array=False):
    # TODO: rewrite the description
    

    nls = len(ls)
    
    mr = cp if use_cp else np # choose matrix routine
        
    assert isinstance(money,tuple)
    assert len(money)==3
    
    
    if use_cp:
        money = tuple((cp.asarray(x) for x in money))
        umult = cp.asarray(umult)
    
    asset_income, wf, wm = money
    
    nexo = wf.size
    na = asset_income.size
    
    
    wf = wf.reshape((1,wf.size))
    wm = wm.reshape((1,wm.size))
    asset_income = asset_income.reshape((asset_income.size,1))
    money = wf + wm + asset_income
        
    
        
    if isinstance(EV,tuple):
        assert len(EV) == 3
        (ind,p,EVin) = (cp.asarray(x) if use_cp else x for x in EV)
        EV_by_l = get_EVM(ind,p,EVin,use_cp)
    
    
    if use_cp: # it is ok to use cp.asarray twice, it does not copy
        money,sgrid,EV_by_l = (cp.asarray(x) for x in (money,sgrid,EV_by_l))
    
    
    ntheta = EV_by_l.shape[-2]
    
    assert money.ndim < EV_by_l.ndim
    assert (EV_by_l.ndim - money.ndim == 2), 'Shape mismatch?'
    shp = money.shape + (ntheta,) # shape of the result
    
    V, c, s = mr.empty(shp,mr.float32), mr.empty(shp,mr.float32), mr.empty(shp,mr.float32)    
    
    oms = 1-sigma
    
    def u(c): return (c**(oms))/(oms)   
    
    ns = sgrid.size
    
    # this will use a weird fact that -2*(1,) = () (empty tuple)
    
    s_size = (1,ns,1)
    
    s_expanded = sgrid.reshape(s_size)
    
    
    tal = cp_take_along_axis if use_cp else np.take_along_axis
    
    
    i_opt_arr = mr.empty((na,nexo,ntheta,nls),dtype=mr.int32)
    c_opt_arr = mr.empty((na,nexo,ntheta,nls),dtype=mr.float32)
    s_opt_arr = mr.empty((na,nexo,ntheta,nls),dtype=mr.float32)
    V_opt_arr = mr.empty((na,nexo,ntheta,nls),dtype=mr.float32)
    
    for i, (lval, uval) in enumerate(zip(ls,us)):
        
        EV_here = EV_by_l[...,i]
        
        money_left = money - (1-lval)*wf.reshape((1,nexo))
        
        if (not use_cp) and (not compare) and (not force_array):
            V, c, s = np.empty((3,na,nexo,ntheta),dtype=np.float32)
            i_opt = -np.ones((na,nexo,ntheta),dtype=np.int16)
                 
            v_opt_couple_local(money_left,sgrid,umult,EV_here,sigma,beta,uval+ushift,V,i_opt,c,s)
        
        else:
            
            if compare:
                V_comp, c_comp, s_comp,  i_opt_comp = (x.copy() for x in (V, c, s, i_opt))
            
            c_mat = mr.expand_dims(money_left,1)  - s_expanded
            u_mat = mr.full(c_mat.shape,-mr.inf)
            u_mat[c_mat>0] = u(c_mat[c_mat>0])
            
            c_mat[c_mat<0] = -1.0
            
            # u_mat shape is (na,ns,nexo)
            u_mat = u_mat.reshape(u_mat.shape + (1,))   # adds dimension for theta 
            u_mat_theta = u_mat*umult.reshape((1,1,1,umult.size)) # we added s dimension :(
            u_mat_total = u_mat_theta + uval + ushift
            
            
            # u_mat_total shape is (na,ns,nexo,ntheta)
            # EV shape is (ns,nexo,ntheta)
            V_arr = u_mat_total + beta*mr.expand_dims(EV_here,0) # adds dimension for current a
            # V_arr shape is (na,ns,nexo,ntheta)
            i_opt = V_arr.argmax(axis=1) # (na,nexo,ntheta)
            
            i_opt_ed = mr.expand_dims(i_opt,1)
            
            
            s = sgrid[i_opt]
            c = tal(mr.expand_dims(c_mat,3),i_opt_ed,1).squeeze(axis=1) #money.reshape( (money.shape+(1,)) ) - s
            V = tal(u_mat_total,i_opt_ed,1).squeeze(axis=1) + beta*tal(EV_here,i_opt,0) # squeeze
        
        
            if compare:
                def maxdiff(x,y): return np.max(np.abs(x-y))
                print('Maximum difference in V is {}'.format(maxdiff(V,V_comp)))
                print('Maximum difference in c is {}'.format(maxdiff(c,c_comp)))
                print('Maximum difference in s is {}'.format(maxdiff(s,s_comp)))
                print('Maximum difference in i_opt is {}'.format(maxdiff(i_opt,i_opt_comp)))
                
                
        i_opt_arr[...,i] = i_opt
        c_opt_arr[...,i] = c
        s_opt_arr[...,i] = s
        V_opt_arr[...,i] = V
        
        
        
        
    i_ls = mr.expand_dims(V_opt_arr.argmax(axis=3),3)
    
    V = tal(V_opt_arr,i_ls,axis=3).squeeze(axis=3)
    c = tal(c_opt_arr,i_ls,axis=3).squeeze(axis=3)
    s = tal(s_opt_arr,i_ls,axis=3).squeeze(axis=3)
    i_opt = tal(i_opt_arr,i_ls,axis=3).squeeze(axis=3)
        
    i_ls = i_ls.squeeze(axis=3)
        
    
    if use_cp:
        ret = lambda x : cp.asnumpy(x)
    else:
        ret = lambda x : x
        
    V_all = V_opt_arr
    
    
    return ret(V), ret(c), ret(s), ret(i_opt), ret(i_ls), ret(V_all).astype(np.float32)






def v_optimize_single(money,sgrid,EV,sigma,beta,ushift,use_cp=ucp,return_ind=False):
    # this is the optimizer for value functions
    # 1. It can use cuda arrays (cupy) if use_cp=True
    # 2. It can accept few shapes of money array and EV
    # 3. If money array and EV array are of different shapes money array is
    # broadcasted (this is for the case where EV varies with theta and money
    # are independent on theta)
    # 4. money can be tuple of two elements (assets income and labor income),
    # this form is more compact. In this case the array if formed by adding
    # them (outter sum). 
    # 5. EV can also be tuple of three elements, that are inputs to get_EVM
    # in this case get_EVM is ran internally, this may help with large arrays
    # so we transition less things to GPU
    # Shape of the result is (money.shape[0],EV.shape[1:])
    
    # TBD: file opt_test.py has jit-able version of these functions .
    # So far they are slower than this but they might be improved
    
    
    
    mr = cp if use_cp else np # choose matrix routine
    
    if isinstance(money,tuple):
        assert len(money) == 2
        
        asset_income = money[0].reshape((money[0].size,1))
        labor_income = money[1].reshape((1,money[1].size))
        if use_cp: asset_income, labor_income = cp.asarray(asset_income), cp.asarray(labor_income)
        money = asset_income + labor_income # broadcasting 
        
    if isinstance(EV,tuple):
        assert len(EV) == 3
        (ind,p,EVin) = (cp.asarray(x) if use_cp else x for x in EV)
        EV = get_EVM(ind,p,EVin,use_cp)
    
    
    if use_cp: # it is ok to use cp.asarray twice, it does not copy
        money,sgrid,EV = (cp.asarray(x) for x in (money,sgrid,EV))
    
    
    assert money.ndim == EV.ndim
        
    assert (money.ndim == EV.ndim), 'Shape mismatch?'
    shp = money.shape
    
    V, c, s = mr.empty(shp,mr.float32), mr.empty(shp,mr.float32), mr.empty(shp,mr.float32)    
    
    oms = 1-sigma
    
    def u(c): return (c**(oms))/(oms)   
    
    ns = sgrid.size
    
    # this will use a weird fact that -2*(1,) = () (empty tuple)
    
    s_size = ((1,ns,1))
    
    s_expanded = sgrid.reshape(s_size)
    
    c_mat = mr.expand_dims(money,1) - s_expanded 
    u_mat = mr.full(c_mat.shape,-mr.inf)
    u_mat[c_mat>0] = u(c_mat[c_mat>0]) + ushift
     
    
    V_arr = u_mat + beta*mr.expand_dims(EV,0)
   
    i_opt = V_arr.argmax(axis=1)
    
    s = sgrid[i_opt]
    
    c = money - s
    
    tal = cp_take_along_axis if use_cp else np.take_along_axis
    
    V = u(c) + beta*tal(EV,i_opt,0)
   
    
    if use_cp:
        ret = lambda x : cp.asnumpy(x)
    else:
        ret = lambda x : x
        
    
    
    if not return_ind:
        return ret(V), ret(c), ret(s)
    else:
        return ret(V), ret(c), ret(s), ret(i_opt)


from numba import prange
@jit(nopython=True,parallel=True)
def v_opt_couple_local(money,sgrid,u_multupliers,EV,sigma,beta,uadd,V_opt,i_opt,c_opt,s_opt):
    # this is a looped version of the optimizer
    # the last two things are outputs
    
    na, nexo, ntheta = money.shape[0], money.shape[1], u_multupliers.size
    
    ns = sgrid.size
    
    assert money.shape == (na,nexo)
    assert V_opt.shape == (na,nexo,ntheta) == i_opt.shape
    assert EV.shape == (ns,nexo,ntheta)
    
    #uout = np.full((na,ns,nexo,ntheta),-np.inf)
    #Vout = np.full((na,ns,nexo,ntheta),-np.inf)
    #cout = np.full((na,ns,nexo),-1.0)
    
    def ufun(x):
        return (x**(1-sigma))/(1-sigma)
    
    for ind_a in prange(na):
        for ind_exo in prange(nexo):
            # finds index of maximum savings
            money_i = money[ind_a,ind_exo]
            if money_i > sgrid[ns-1]: # if can save the max amount
                ind_s = ns-1 
            else:
                # if not look for the highest feasible savings
                for ind_s in range(ns-1):
                    if money_i <= sgrid[ind_s+1]: break
                assert money_i > sgrid[ind_s]
                
            
            
            # finds utility values of all feasible savings
            # should I make a loop instead?
            u_of_s = ufun( money_i - sgrid[0:ind_s+1] )
            
            # each value of theta corresponds to different EV
            # so do the maximum search for each theta using these utilities
            
            for ind_theta in prange(ntheta):
                mult = u_multupliers[ind_theta]
                EV_of_s = EV[0:(ind_s+1),ind_exo,ind_theta]                
                u_adjusted = mult*u_of_s + uadd
                V_of_s = u_adjusted + beta*EV_of_s
                
                
                io = V_of_s.argmax()
                i_opt[ind_a,ind_exo,ind_theta] = io
                V_opt[ind_a,ind_exo,ind_theta] = V_of_s[io] 
                c_opt[ind_a,ind_exo,ind_theta] = money_i - sgrid[io]
                s_opt[ind_a,ind_exo,ind_theta] = sgrid[io]
    
    
                
            