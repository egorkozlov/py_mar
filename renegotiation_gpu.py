#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 19:56:40 2020

@author: egorkozlov
"""

from numba import cuda, f4, f8
import numpy as np

use_f32 = False

if use_f32:
    gpu_type = f4
    cpu_type = np.float32
else:
    gpu_type = f8
    cpu_type = np.float64


def v_ren_gpu_oneopt(v_y_ni, vf_y_ni, vm_y_ni, vf_n_ni, vm_n_ni, itht, wntht, thtgrid, 
                          rescale = True):
    
                    
    na, ne, nt_coarse = v_y_ni.shape
    nt = thtgrid.size
    assert rescale, 'no rescale is not implemented'
    
    assert nt < 500 
    
    
    
    v_out = np.empty((na,ne,nt),dtype=cpu_type)
    vm_out = np.empty((na,ne,nt),dtype=cpu_type)
    vf_out = np.empty((na,ne,nt),dtype=cpu_type)
    itheta_out = np.empty((na,ne,nt),dtype=np.int16)
    
    thtgrid = cuda.to_device(thtgrid)
    
    for ia in range(na):
    
        threadsperblock = (1, nt)
            
        b_exo = ne
        b_theta = 1
        
        blockspergrid = (b_exo, b_theta)
        
        v_yi, vf_yi, vm_yi = [cuda.to_device(
                                np.ascontiguousarray(x[ia,:,:].copy())
                                            ) for x in (v_y_ni, vf_y_ni, vm_y_ni)]
        
        vf_ni, vm_ni = [cuda.to_device(
                                        np.ascontiguousarray(x[ia,:,:].copy())
                                      ) for x in (vf_n_ni,vm_n_ni)]
        
        
        v_outi = cuda.device_array((ne,nt),dtype=cpu_type)
        vm_outi = cuda.device_array((ne,nt),dtype=cpu_type)
        vf_outi = cuda.device_array((ne,nt),dtype=cpu_type)
        itheta_outi = cuda.device_array((ne,nt),dtype=np.int16)
        
                         
        
        cuda_ker_one_opt[blockspergrid, threadsperblock](v_yi, vf_yi, vm_yi, vf_ni, vm_ni, 
                                        itht, wntht, thtgrid,  
                                        v_outi, vm_outi, vf_outi, itheta_outi)
    
    
        v_out[ia,:,:] = v_outi
        vm_out[ia,:,:] = vm_outi
        vf_out[ia,:,:] = vf_outi
        itheta_out[ia,:,:] = itheta_outi
    
    return v_out, vf_out, vm_out, itheta_out
    



@cuda.jit   
def cuda_ker_one_opt(v_y_ni, vf_y_ni, vm_y_ni, vf_n_ni, vm_n_ni, itht, wntht, thtgrid, v_out, vm_out, vf_out, itheta_out):
    # this assumes block is for the same a and theta
    ie, it = cuda.grid(2)
    
    v_in_store  = cuda.shared.array((500,),gpu_type)
    vf_in_store = cuda.shared.array((500,),gpu_type)
    vm_in_store = cuda.shared.array((500,),gpu_type)
    
    vf_no_store = cuda.shared.array((500,),gpu_type)
    vm_no_store = cuda.shared.array((500,),gpu_type)
    
    
    
    ne = v_y_ni.shape[0]
    nt_crude = v_y_ni.shape[1]
    nt = thtgrid.size
    
    
    f1 = gpu_type(1.0)
    
    if ie < ne and it < nt:
        
        it_int = itht[it]
        for ittc in range(nt_crude):
            if ittc==it_int: break
        
        
        ittp = ittc + 1
        wttp = wntht[it]
        wttc = f1 - wttp
        
        
        v_in_store[it]  = wttc*v_y_ni[ie,ittc]  + wttp*v_y_ni[ie,ittp]
        vf_in_store[it] = wttc*vf_y_ni[ie,ittc] + wttp*vf_y_ni[ie,ittp]
        vm_in_store[it] = wttc*vm_y_ni[ie,ittc] + wttp*vm_y_ni[ie,ittp]
        
        vf_no_store[it] = wttc*vf_n_ni[ie,ittc] + wttp*vf_n_ni[ie,ittp]
        vm_no_store[it] = wttc*vm_n_ni[ie,ittc] + wttp*vm_n_ni[ie,ittp]
        
        
        v_out[ie,it] = v_in_store[it]
        vf_out[ie,it] = vf_in_store[it] 
        vm_out[ie,it] = vm_in_store[it] 
        
        cuda.syncthreads()
        
        vf_no = vf_no_store[it]
        vm_no = vm_no_store[it]
        
        
        
        
        
        # fill default values
        
        
        
        if vf_in_store[it] >= vf_no and vm_in_store[it] >= vm_no:
        #if vf_y[ie,it] >= vf_no and vm_y[ie,it] >= vm_no:
            itheta_out[ie,it] = it
            return
        
        if vf_in_store[it] < vf_no and vm_in_store[it] < vm_no:
        #if vf_y[ie,it] < vf_no and vm_y[ie,it] < vm_no:
            itheta_out[ie,it] = -1
            tht = thtgrid[it]
            v_out[ie,it] = tht*vf_no + (1-tht)*vm_no
            vf_out[ie,it] = vf_no
            vm_out[ie,it] = vm_no
            return
        
        
       
        it_ren = -1
        
        found_increase = False
        found_decrease = False
         
         
        for it_inc in range(it+1,nt):
            if (vf_in_store[it_inc] >= vf_no and vm_in_store[it_inc] >= vm_no):
                found_increase = True
                break
         
        for it_dec in range(it-1,-1,-1):
            if (vf_in_store[it_dec] >= vf_no and vm_in_store[it_dec] >= vm_no):
                found_decrease = True
                break
#            
#        
        if found_increase and found_decrease:
            dist_increase = it_inc - it
            dist_decrease = it - it_dec
#            
            if dist_increase != dist_decrease:
                it_ren = it_inc if dist_increase < dist_decrease else it_dec
            else:
                # tie breaker
                # numba-cuda does not do abs so we do these dumb things
                dist_mid_inc = it_inc - (nt/2)                
                if dist_mid_inc < 0: dist_mid_inc = -dist_mid_inc
                dist_mid_dec = it_dec - (nt/2)
                if dist_mid_dec < 0: dist_mid_dec = -dist_mid_dec
                it_ren = it_inc if dist_mid_inc < dist_mid_dec else it_dec
            
        elif found_increase and not found_decrease:
            it_ren = it_inc
        elif found_decrease and not found_increase:
            it_ren = it_dec
        else:
            it_ren = -1 # check this!
             
         # finally fill the values    
             
        if it_ren == -1:
            tht = thtgrid[it]
            v_out[ie,it] = tht*vf_no + (1-tht)*vm_no
            vf_out[ie,it] = vf_no
            vm_out[ie,it] = vm_no
            itheta_out[ie,it] = -1
        else:
             
            # rescaling
            tht_old = thtgrid[it]
            tht_new = thtgrid[it_ren]
            factor = (1-tht_old)/(1-tht_new) if tht_old < tht_new else tht_old/tht_new
             
            v_out[ie,it] = factor*v_in_store[it_ren]
            vf_out[ie,it] = vf_in_store[it_ren]
            vm_out[ie,it] = vm_in_store[it_ren]
            itheta_out[ie,it] = it_ren





def v_ren_gpu_twoopt(v_y_ni0, v_y_ni1, vf_y_ni0, vf_y_ni1, vm_y_ni0, vm_y_ni1, vf_n_ni, vm_n_ni, itht, wntht, thtgrid, 
                          rescale = True):
    
                    
    na, ne, nt_coarse = v_y_ni0.shape
    nt = thtgrid.size
    assert rescale, 'no rescale is not implemented'
    
    assert nt < 500 
    
    
    
    
    v_out = np.empty((na,ne,nt),dtype=cpu_type)
    vm_out = np.empty((na,ne,nt),dtype=cpu_type)
    vf_out = np.empty((na,ne,nt),dtype=cpu_type)
    itheta_out = np.empty((na,ne,nt),dtype=np.int16)
    switch_out = np.empty((na,ne,nt),dtype=np.bool_)
    
    thtgrid = cuda.to_device(thtgrid)
    
    for ia in range(na):
    
        threadsperblock = (1, nt)
            
        b_exo = ne
        b_theta = 1
        
        blockspergrid = (b_exo, b_theta)
        
        v_yi0, vf_yi0, vm_yi0 = [cuda.to_device(
                                np.ascontiguousarray(x[ia,:,:].copy())
                                            ) for x in (v_y_ni0, vf_y_ni0, vm_y_ni0)]
    
        v_yi1, vf_yi1, vm_yi1 = [cuda.to_device(
                                np.ascontiguousarray(x[ia,:,:].copy())
                                            ) for x in (v_y_ni1, vf_y_ni1, vm_y_ni1)]
        
        vf_ni, vm_ni = [cuda.to_device(
                                        np.ascontiguousarray(x[ia,:,:].copy())
                                      ) for x in (vf_n_ni,vm_n_ni)]
        
        
        v_outi = cuda.device_array((ne,nt),dtype=cpu_type)
        vm_outi = cuda.device_array((ne,nt),dtype=cpu_type)
        vf_outi = cuda.device_array((ne,nt),dtype=cpu_type)
        itheta_outi = cuda.device_array((ne,nt),dtype=np.int16)
        switch_outi = cuda.device_array((ne,nt),dtype=np.bool_)
        
                         
        
        cuda_ker_two_opt[blockspergrid, threadsperblock](v_yi0, v_yi1, vf_yi0, vf_yi1, vm_yi0, vm_yi1, vf_ni, vm_ni, 
                                        itht, wntht, thtgrid,  
                                        v_outi, vm_outi, vf_outi, itheta_outi, switch_outi)
    
    
        v_out[ia,:,:] = v_outi
        vm_out[ia,:,:] = vm_outi
        vf_out[ia,:,:] = vf_outi
        itheta_out[ia,:,:] = itheta_outi
        switch_out[ia,:,:] = switch_outi
    
    return v_out, vf_out, vm_out, itheta_out, switch_out
    


            
            

@cuda.jit   
def cuda_ker_two_opt(v_y_ni0, v_y_ni1, vf_y_ni0, vf_y_ni1, vm_y_ni0, vm_y_ni1, vf_n_ni, vm_n_ni, itht, wntht, thtgrid, v_out, vm_out, vf_out, itheta_out, switch_out):
    # this assumes block is for the same a and theta
    ie, it = cuda.grid(2)
    
    v_in_store  = cuda.shared.array((500,),gpu_type)
    vf_in_store = cuda.shared.array((500,),gpu_type)
    vm_in_store = cuda.shared.array((500,),gpu_type)
    
    vf_no_store = cuda.shared.array((500,),gpu_type)
    vm_no_store = cuda.shared.array((500,),gpu_type)
    
    
    
    ne = v_y_ni0.shape[0]
    nt_crude = v_y_ni0.shape[1]
    nt = thtgrid.size
    
    
    f1 = f8(1.0)
    
    if ie < ne and it < nt:
        
        it_int = itht[it]
        for ittc in range(nt_crude):
            if ittc==it_int: break
        
        
        ittp = ittc + 1
        wttp = wntht[it]
        wttc = f1 - wttp
        
        def wsum(x):
            return wttc*x[ie,ittc] + wttp*x[ie,ittp]
        
        vy_0 = wsum(v_y_ni0)
        vy_1 = wsum(v_y_ni1)
        
        pick1 = (vy_1 > vy_0)
        #print(vy_1-vy_0)
        
        switch_out[ie,it] = pick1
        
        if pick1:
            v_in_store[it]  = vy_1
            vf_in_store[it] = wsum(vf_y_ni1)
            vm_in_store[it] = wsum(vm_y_ni1)
        else:
            v_in_store[it]  = vy_0
            vf_in_store[it] = wsum(vf_y_ni0)
            vm_in_store[it] = wsum(vm_y_ni0)
        
        
        vf_no_store[it] = wttc*vf_n_ni[ie,ittc] + wttp*vf_n_ni[ie,ittp]
        vm_no_store[it] = wttc*vm_n_ni[ie,ittc] + wttp*vm_n_ni[ie,ittp]
        
        
        
        v_out[ie,it] = v_in_store[it]
        vf_out[ie,it] = vf_in_store[it] 
        vm_out[ie,it] = vm_in_store[it] 
        
        cuda.syncthreads()
        
        vf_no = vf_no_store[it]
        vm_no = vm_no_store[it]
        
        
        
        if vf_in_store[it] >= vf_no and vm_in_store[it] >= vm_no:
        #if vf_y[ie,it] >= vf_no and vm_y[ie,it] >= vm_no:
            itheta_out[ie,it] = it
            return
        
        if vf_in_store[it] < vf_no and vm_in_store[it] < vm_no:
        #if vf_y[ie,it] < vf_no and vm_y[ie,it] < vm_no:
            itheta_out[ie,it] = -1
            tht = thtgrid[it]
            v_out[ie,it] = tht*vf_no + (1-tht)*vm_no
            vf_out[ie,it] = vf_no
            vm_out[ie,it] = vm_no
            return
        
        
       
        it_ren = -1
        
        found_increase = False
        found_decrease = False
         
         
        for it_inc in range(it+1,nt):
            if (vf_in_store[it_inc] >= vf_no and vm_in_store[it_inc] >= vm_no):
                found_increase = True
                break
         
        for it_dec in range(it-1,-1,-1):
            if (vf_in_store[it_dec] >= vf_no and vm_in_store[it_dec] >= vm_no):
                found_decrease = True
                break
#            
#        
        if found_increase and found_decrease:
            dist_increase = it_inc - it
            dist_decrease = it - it_dec
#            
            if dist_increase != dist_decrease:
                it_ren = it_inc if dist_increase < dist_decrease else it_dec
            else:
                # tie breaker
                # numba-cuda does not do abs so we do these dumb things
                dist_mid_inc = it_inc - (nt/2)                
                if dist_mid_inc < 0: dist_mid_inc = -dist_mid_inc
                dist_mid_dec = it_dec - (nt/2)
                if dist_mid_dec < 0: dist_mid_dec = -dist_mid_dec
                it_ren = it_inc if dist_mid_inc < dist_mid_dec else it_dec
            
        elif found_increase and not found_decrease:
            it_ren = it_inc
        elif found_decrease and not found_increase:
            it_ren = it_dec
        else:
            it_ren = -1 # check this!
             
         # finally fill the values    
             
        if it_ren == -1:
            tht = thtgrid[it]
            v_out[ie,it] = tht*vf_no + (1-tht)*vm_no
            vf_out[ie,it] = vf_no
            vm_out[ie,it] = vm_no
            itheta_out[ie,it] = -1
        else:
             
            # rescaling
            tht_old = thtgrid[it]
            tht_new = thtgrid[it_ren]
            factor = (1-tht_old)/(1-tht_new) if tht_old < tht_new else tht_old/tht_new
             
            v_out[ie,it] = factor*v_in_store[it_ren]
            vf_out[ie,it] = vf_in_store[it_ren]
            vm_out[ie,it] = vm_in_store[it_ren]
            itheta_out[ie,it] = it_ren