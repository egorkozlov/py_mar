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


if system() != 'Darwin' and system() != 'Windows':
#if system() != 'Darwin':
    import cupy as cp
    ugpu = True
else:
    ugpu = False


#@jit('float64[:](float64[:], int64, float64, float64)',nopython=True)



def get_EVM(ind,wthis,EVin,use_gpu=False,dtype=np.float32):
    # this essenitally doubles VecOnGrid.apply method
    # but we cannot deprecate it unless we are sure VecOnGrid works on GPU 
    # correctly
    
    mr = cp if use_gpu else np
    
    
    ev_aux_shape  = EVin.shape[1:]
    shap = (ind.size,) + ev_aux_shape
    
    
    if dtype == np.float32:
        EVout = mr.empty(shap,mr.float32)
    else:
        if use_gpu:
            assert dtype==np.float64
            EVout = mr.empty(shap,mr.float64)
        else:
            EVout = np.empty(shap,dtype)
    
    pb = wthis.reshape(((wthis.size,)+(1,)*len(ev_aux_shape)))
    EVout[:] = pb*EVin[ind,...] + (1-pb)*EVin[ind+1,...]    
    
    return EVout



def v_optimize_couple(money_in,sgrid,EV,mgrid,utilint,xint,ls,beta,ushift,use_gpu=ugpu,dtype=np.float32):
    # This optimizer avoids creating big arrays and uses parallel-CPU on 
    # machines without NUMBA-CUDA codes otherwise
    

    nls = len(ls)
    
        
    assert isinstance(money_in,tuple)
    assert len(money_in)==3
    
    asset_income, wf, wm = money_in
    
    
    nexo = wf.size
    na = asset_income.size
    
    wf = wf.reshape((1,wf.size))
    wm = wm.reshape((1,wm.size))
    asset_income = asset_income.reshape((asset_income.size,1))
    money = wf + wm + asset_income
        
    
    if isinstance(EV,tuple):
        assert len(EV) == 3
        ind,p,EVin = EV
        EV_by_l = get_EVM(ind,p,EVin,use_gpu=False,dtype=dtype) 
    
    ntheta = EV_by_l.shape[-2]
    
    tal = np.take_along_axis
    
    i_opt_arr = np.empty((na,nexo,ntheta,nls),dtype=np.int16)
    c_opt_arr = np.empty((na,nexo,ntheta,nls),dtype=dtype)
    x_opt_arr = np.empty((na,nexo,ntheta,nls),dtype=dtype)
    s_opt_arr = np.empty((na,nexo,ntheta,nls),dtype=dtype)
    V_opt_arr = np.empty((na,nexo,ntheta,nls),dtype=dtype)
    
    for i, lval in enumerate(ls):
        
        EV_here = EV_by_l[...,i]
        util = utilint[...,i]
        xvals = xint[...,i]
        
        money_left = money - (1-lval)*wf.reshape((1,nexo))
        
        if not use_gpu:
            
            # preallocation helps a bit here          
            V, c, x, s = np.empty((4,na,nexo,ntheta),dtype=dtype)
            i_opt = -np.ones((na,nexo,ntheta),dtype=np.int16)                 
            v_couple_local_intu(money_left,sgrid,EV_here,mgrid,util,xvals,beta,ushift,V,i_opt,c,x,s)
            
            V2, i_opt2, c2, x2, s2 = v_couple_gpu(money_left,sgrid,EV_here,mgrid,util,xvals,beta,ushift)

            assert np.allclose(V2,V)
            assert np.all(i_opt==i_opt2)
            
        else:
            #print((EV_here.shape,(na,nexo,ntheta,nls)))
            V, i_opt, c, x, s = v_couple_gpu(money_left,sgrid,EV_here,mgrid,util,xvals,beta,ushift)


        
                
        i_opt_arr[...,i] = i_opt
        c_opt_arr[...,i] = c
        x_opt_arr[...,i] = x
        s_opt_arr[...,i] = s
        V_opt_arr[...,i] = V # discard this V later as it is not very precise
        
        
        
    i_ls = np.expand_dims(V_opt_arr.argmax(axis=3),3)
    
    V = tal(V_opt_arr,i_ls,axis=3).squeeze(axis=3)
    c = tal(c_opt_arr,i_ls,axis=3).squeeze(axis=3)
    s = tal(s_opt_arr,i_ls,axis=3).squeeze(axis=3)
    i_opt = tal(i_opt_arr,i_ls,axis=3).squeeze(axis=3)
        
    i_ls = i_ls.squeeze(axis=3)
        
    
    ret = lambda q : q
        
    V_all = V_opt_arr
    
    
    
    return ret(V), ret(c), ret(x), ret(s), ret(i_opt), ret(i_ls), ret(V_all).astype(dtype)







def v_optimize_single(money,sgrid,EV,sigma,beta,ushift,use_gpu=False,return_ind=False,dtype=np.float32):
    # this is the optimizer for value functions
    # 1. It can use cuda arrays (cupy) if use_gpu=True
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
    
    
    
    mr = cp if use_gpu else np # choose matrix routine
    
    if use_gpu:
        if dtype==np.float32:
            dtype_here = mr.float32 
        elif dtype==np.float64:
            dtype_here = mr.float64
        else:
            raise(Exception('unsupported type...'))
    else:
        dtype_here = dtype
    
    if isinstance(money,tuple):
        assert len(money) == 2
        
        asset_income = money[0].reshape((money[0].size,1))
        labor_income = money[1].reshape((1,money[1].size))
        if use_gpu: asset_income, labor_income = cp.asarray(asset_income), cp.asarray(labor_income)
        money = asset_income + labor_income # broadcasting 
        
    if isinstance(EV,tuple):
        assert len(EV) == 3
        (ind,p,EVin) = (cp.asarray(x) if use_gpu else x for x in EV)
        EV = get_EVM(ind,p,EVin,use_gpu,dtype=dtype)
    
    
    if use_gpu: # it is ok to use cp.asarray twice, it does not copy
        money,sgrid,EV = (cp.asarray(x) for x in (money,sgrid,EV))
    
    
    assert money.ndim == EV.ndim
        
    assert (money.ndim == EV.ndim), 'Shape mismatch?'
    shp = money.shape
    
    V, c, s = mr.empty(shp,dtype_here), mr.empty(shp,dtype_here), mr.empty(shp,dtype_here)    
    
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
    
    tal = cp_take_along_axis if use_gpu else np.take_along_axis
    
    V = u(c) + beta*tal(EV,i_opt,0)
   
    
    if use_gpu:
        ret = lambda x : cp.asnumpy(x)
    else:
        ret = lambda x : x
        
    
    
    if not return_ind:
        return ret(V), ret(c), ret(s)
    else:
        return ret(V), ret(c), ret(s), ret(i_opt)


from numba import prange                
@jit(nopython=True)#,parallel=True)
def v_couple_local_intu(money,sgrid,EV,mgrid,u_on_mgrid,x_on_mgrid,beta,uadd,V_opt,i_opt,c_opt,x_opt,s_opt):
    # this is a looped version of the optimizer
    # the last two things are outputs
    
    na, nexo, ntheta = money.shape[0], money.shape[1], EV.shape[2]
    
    ns = sgrid.size
    
    nm = mgrid.size
    
    assert money.shape == (na,nexo)
    assert V_opt.shape == (na,nexo,ntheta) == i_opt.shape
    assert EV.shape == (ns,nexo,ntheta)
    
    
    
    coh_min = mgrid[0] 
    
    for ind_a in prange(na):
        for ind_exo in prange(nexo):
            # finds index of maximum savings
            money_i = money[ind_a,ind_exo]
            money_minus_coh = money_i - coh_min
            
            ind_s = np.minimum( np.searchsorted(sgrid,money_minus_coh)-1,ns-1)
            
            i_m = np.minimum( np.searchsorted(mgrid,money_i)-1,nm-1)
            
            i_m_all = np.zeros((ind_s+1,),np.int16)
            
            i_m_all[0] = i_m
            
            for i_cand in range(1,ind_s+1):
                m_after_s = money_i - sgrid[i_cand]
                while mgrid[i_m] > m_after_s:
                    i_m -= 1                    
                assert i_m >= 0
                i_m_all[i_cand] = i_m
                
            
            for ind_theta in range(ntheta):
                
                ugrid_ce = u_on_mgrid[:,ind_theta]
                bEVval = beta*EV[:,ind_exo,ind_theta]
                
                io = 0
                Vo = -1e20
                
                
                for i_cand in range(ind_s+1):
                    
                    i_m = i_m_all[i_cand]
                    u_cand = ugrid_ce[i_m]
                    
                    V_cand = u_cand + bEVval[i_cand]
                    
                    if i_cand == 0 or V_cand > Vo:
                        io, Vo = i_cand, V_cand                     
                    
                
                i_opt[ind_a,ind_exo,ind_theta] = io
                V_opt[ind_a,ind_exo,ind_theta] = Vo + uadd# NB: this V is imprecise
                # you can recover V from optimal savings & consumption later
                x = x_on_mgrid[i_m_all[io],ind_theta]
                x_opt[ind_a,ind_exo,ind_theta] = x
                c_opt[ind_a,ind_exo,ind_theta] = money_i - x - sgrid[io]
                s_opt[ind_a,ind_exo,ind_theta] = sgrid[io]        
                assert Vo > -1e20
            
    

from math import ceil
from numba import cuda

def v_couple_gpu(money,sgrid,EV,mgrid,util,xvals,beta,uadd,use_kernel_pool=False):
    
    
    na, nexo, ntheta = money.shape[0], money.shape[1], EV.shape[2]
    
    ns = sgrid.size
    
    
    #assert ns < 2000, 'Please alter the array size in cuda_ker_pool'
    
    assert money.shape == (na,nexo)
    assert EV.shape == (ns,nexo,ntheta)
    
    V_opt_g = cuda.device_array((na,nexo,ntheta),dtype=EV.dtype)  
    x_opt_g = cuda.device_array((na,nexo,ntheta),dtype=EV.dtype)  
    i_opt_g = cuda.device_array((na,nexo,ntheta),dtype=np.int16)
    
    bEV = beta*EV
    
    
    money_g, sgrid_g, bEV_g, mgrid_g, util_g, xvals_g = (cuda.to_device(np.ascontiguousarray(x)) for x in (money, sgrid, bEV, mgrid, util, xvals))
    
    
    threadsperblock = (8, 16, 8)
    # this is a tunning parameter. 8*16*8=1024 is the number of threads per 
    # block. This is GPU specific, on Quests's GPU the maximum is 1024, on 
    # different machines it can be lower
    
    b_a = int(ceil(na / threadsperblock[0]))
    b_exo = int(ceil(nexo / threadsperblock[1]))
    b_theta = int(ceil(ntheta / threadsperblock[2]))
    blockspergrid = (b_a, b_exo, b_theta)
    
    cuda_ker[blockspergrid, threadsperblock](money_g, sgrid_g, bEV_g, mgrid_g, util_g, xvals_g, uadd,
                                                V_opt_g,i_opt_g,x_opt_g)
    
    V_opt, i_opt, x_opt = (x.copy_to_host() for x in (V_opt_g,i_opt_g,x_opt_g))
    
    s_opt = sgrid[i_opt]
    c_opt = money[:,:,None] - x_opt - s_opt
    
    return V_opt,i_opt,c_opt,x_opt,s_opt


#from numba import f4

@cuda.jit
def cuda_ker(money_g, sgrid_g, bEV_g, mgrid_g, util_g, xvals_g, uadd, V_opt_g, i_opt_g, x_opt_g):
    ind_a, ind_exo, ind_theta = cuda.grid(3)
        
    na = money_g.shape[0]
    nexo = money_g.shape[1]
    ntheta = bEV_g.shape[2]
    ns = sgrid_g.size
    nm = mgrid_g.size
    
    if (ind_a < na) and (ind_exo < nexo) and (ind_theta < ntheta):
        # finds index of maximum savings
        money_i = money_g[ind_a,ind_exo]
        bEV_of_s = bEV_g[:,ind_exo,ind_theta]  
        util_here = util_g[:,ind_theta]
        xvals_here = xvals_g[:,ind_theta]
        
        mmin = mgrid_g[0]
        
            
        for im in range(nm-1,-1,-1):
            if mgrid_g[im] < money_i: break
        
        iobest = 0
        vbest = util_here[im] + bEV_of_s[0]
        xbest = xvals_here[im]
        
        for io in range(1,ns):
            money_left = money_i - sgrid_g[io]
            if money_left < mmin: break
            
            while mgrid_g[im] > money_left:
                im -= 1
            
            assert im >= 0            
            
            vnow = util_here[im] + bEV_of_s[io] 
            if vnow > vbest:
                iobest = io
                vbest = vnow
                xbest = xvals_here[im] 
            
        i_opt_g[ind_a,ind_exo,ind_theta] = iobest
        x_opt_g[ind_a,ind_exo,ind_theta] = xbest
        V_opt_g[ind_a,ind_exo,ind_theta] = vbest + uadd
        
    
    
