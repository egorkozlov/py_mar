#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects selected routines from opt_test.py

The aim is keep each function readable, and either njit-able or using cuda

Both versions are needed as it is not clear whether we'll actually use cuda

"""

import numpy as np
from numba import jit, prange, cuda, float32
from platform import system
from aux_routines import cp_take_along_axis

if system() != 'Darwin':
    import cupy as cp
    ucp = True
else:
    ucp = False


@jit('float64[:](float64[:], int64, float64, float64)',nopython=True)
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


@jit('Tuple((int32[:],float32[:]))(float64[:],float64[:])',nopython=True)
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
def get_EVM(ind,p,EVin):
    
    ev_aux_shape  = EVin.shape[1:]
    shap = (ind.size,) + ev_aux_shape
    EVout = np.empty(shap,np.float32)
    
    pb = p.reshape(((p.size,)+(1,)*len(ev_aux_shape)))
    EVout[:] = pb*EVin[ind,...] + (1-pb)*EVin[ind+1,...]    
    return EVout


@jit(nopython=True,parallel=True)
def v_optimize_1d_loop(money,sgrid,EV_on_sgrid_1d,sigma,beta):
    # this assumes EV to be one-dimensional
    # this can provide some speedup if the loop is not too long
    
    # it can be two cases: either we need to broadcast u or we do not
    # maybe compute u first and then do broadcasting
    
    
    V, c, s = np.empty(money.shape,np.float64), np.empty(money.shape,np.float64), np.empty(money.shape,np.float64)
    
    ns = sgrid.size
    
    oms = 1-sigma
    
    # does not allow for logs yet
    def u(c): return (c**(oms))/(oms)    
    
    for im in range(money.size): # or prange
        mi = money[im]        
        c_best = mi
        V_best = np.float32(-np.inf)
        s_best = np.float32(0.0)   
        
        for i_sav in range(ns):
            s_current = sgrid[i_sav]
            c_current = mi - s_current
            if c_current <= 0: break                 
            V_new = u(c_current) + beta*EV_on_sgrid_1d[i_sav]
            
            # update if imporved
            if V_new >= V_best:
                c_best = c_current
                V_best = V_new
                s_best = s_current
                    
        c[im] = c_best
        s[im] = s_best
        V[im] = V_best
        
    return V, c, s

def v_optimize(money,sgrid,EVM,sigma,beta,use_cp=ucp):
    
    
    
    
    mr = cp if use_cp else np
    
    if use_cp:
        money,sgrid,EVM = (cp.asarray(x) for x in (money,sgrid,EVM))
    
    
    ntheta = EVM.shape[-1]
    
    if money.ndim < EVM.ndim:
        shp = money.shape + (ntheta,) # shape of the result
        broadcast_um = True
    else:
        assert money.ndim == EVM.ndim, 'Shape missmatch?'
        shp = money.shape
        broadcast_um = False
    
    V, c, s = mr.empty(shp,mr.float32), mr.empty(shp,mr.float32), mr.empty(shp,mr.float32)    
    
    oms = 1-sigma
    
    def u(c): return (c**(oms))/(oms)   
    
    ns = sgrid.size
    
    # this will use a weird fact that -2*(1,) = () (empty tuple)
    
    s_size = ((1,ns) + (EVM.ndim-2)*(1,))
    
    s_expanded = sgrid.reshape(s_size)
    
    c_mat = mr.expand_dims(money,1) - s_expanded 
    u_mat = mr.full(c_mat.shape,-mr.inf)
    u_mat[c_mat>0] = u(c_mat[c_mat>0])
    
    if broadcast_um:
        u_mat = u_mat.reshape(u_mat.shape + (1,))    
    
    V_arr = u_mat + beta*mr.expand_dims(EVM,0)
    
    i_opt = V_arr.argmax(axis=1)
    
    s = sgrid[i_opt]
    
    c = money.reshape( (money.shape+(EVM.ndim-money.ndim)*(1,)) ) - s
    
    tal = cp_take_along_axis if use_cp else np.take_along_axis
    
    V = u(c) + beta*tal(EVM,i_opt,0)
    
    if use_cp:
        ret = lambda x : cp.asnumpy(x)
    else:
        ret = lambda x : x
    
    return ret(V), ret(c), ret(s)

