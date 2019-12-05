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


#@jit('Tuple((int32[:],float32[:]))(float64[:],float64[:])',nopython=True)
def sgrid_on_agrid(sgrid,agrid):
    
    ind = np.empty(sgrid.size,dtype=np.int32)
    p = np.empty(sgrid.size,dtype=np.float32)
    # this uses the fact that both agrid and sgrid are sorted
    # this returns indices and weights that are needed for interpolation of 
    # something that is defined on agrid for values of sgrid
    
    # sgrid[i] = p[i]*agrid[ind[i]] + (1-p[i])*agrid[ind[i]+1]
    
    na = agrid.size
    
    i_next_a = 1
    a_next, a_this = agrid[1], agrid[0]
    
    for i in range(sgrid.size):
        s_this = sgrid[i]
        while a_next < s_this and i_next_a < na-1:
            a_this = agrid[i_next_a]
            i_next_a += 1
            a_next = agrid[i_next_a]
    
        
        #assert s_this >= a_this
        
        ind[i] = i_next_a - 1        
        p[i] = (a_next - s_this)/(a_next - a_this)
        
    return ind, p


# the signature is not valid as EVin may be of variable size, multiple
# signatures are advised if we need providing signatures
#@jit('float32[:](int32[:],float32[:],float32[:])',nopython=True)
def get_EVM(ind,p,EVin,use_cp=False):
    
    mr = cp if use_cp else np
    
    ev_aux_shape  = EVin.shape[1:]
    shap = (ind.size,) + ev_aux_shape
    EVout = mr.empty(shap,mr.float32)
    
    pb = p.reshape(((p.size,)+(1,)*len(ev_aux_shape)))
    EVout[:] = pb*EVin[ind,...] + (1-pb)*EVin[ind+1,...]    
    return EVout


    

def v_optimize_couple(money,sgrid,umult,EV,sigma,beta,ls,us,use_cp=ucp):
    # TODO: rewrite the description
    

    nls = len(ls)
    
    mr = cp if use_cp else np # choose matrix routine
        
    assert isinstance(money,tuple)
    assert len(money)==3
    
    
    if use_cp:
        money = tuple((cp.asarray(x) for x in money))
        umult = cp.asarray(umult)
    
    
    
    wf = money[1]
    wm = money[2]
    asset_income = money[0]
    
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
        
        c_mat = mr.expand_dims(money,1) - s_expanded - (1-lval)*wf.reshape((1,1,nexo))
        u_mat = mr.full(c_mat.shape,-mr.inf)
        u_mat[c_mat>0] = u(c_mat[c_mat>0])
        # u_mat shape is (na,ns,nexo)
        u_mat = u_mat.reshape(u_mat.shape + (1,))   # adds dimension for theta 
        u_mat_theta = u_mat*umult.reshape((1,1,1,umult.size)) # we added s dimension :(
        u_mat_total = u_mat_theta + uval
        
        
        # u_mat_total shape is (na,ns,nexo,ntheta)
        # EV shape is (ns,nexo,ntheta)
        V_arr = u_mat_theta + beta*mr.expand_dims(EV_here,0) # adds dimension for current a
        # V_arr shape is (na,ns,nexo,ntheta)
        i_opt = V_arr.argmax(axis=1) # (na,nexo,ntheta)
        
        i_opt_ed = mr.expand_dims(i_opt,1)
        
        
        s = sgrid[i_opt]
        c = tal(mr.expand_dims(c_mat,3),i_opt_ed,1).squeeze(axis=1) #money.reshape( (money.shape+(1,)) ) - s
        
        V = tal(u_mat_total,i_opt_ed,1).squeeze(axis=1) + beta*tal(EV_here,i_opt,0) # squeeze
        
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
    
    
    return ret(V), ret(c), ret(s), ret(i_opt), ret(i_ls), ret(V_all)




def v_optimize_single(money,sgrid,EV,sigma,beta,use_cp=ucp,return_ind=False):
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
    u_mat[c_mat>0] = u(c_mat[c_mat>0])
     
    
    V_arr = u_mat + beta*mr.expand_dims(EV,0)
   
    i_opt = V_arr.argmax(axis=1)
    
    s = sgrid[i_opt]
    
    c = money - s
    
    tal = cp_take_along_axis if use_cp else np.take_along_axis
    
    V = u(c) + beta*tal(EV,i_opt,0)
    #V2 = V_arr.max(axis=1)
    #print('Maximum difference is {}'.format(np.max(np.abs(V2-V))))
    
    
    
    if use_cp:
        ret = lambda x : cp.asnumpy(x)
    else:
        ret = lambda x : x
        
    
    
    if not return_ind:
        return ret(V), ret(c), ret(s)
    else:
        return ret(V), ret(c), ret(s), ret(i_opt)

    

'''

'''
