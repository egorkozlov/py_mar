#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This tests routine about finding the points where thing hits zero with respect
to the last dimension of the array
"""

import numpy as np

# this is a dumb non-vectorized thing
def zero_hit_loop(vin):
    assert vin.ndim == 1
    k = np.empty((vin.size-1,),dtype=np.float64)
    
    for ik in range(vin.size-1):
        dvi = vin[ik+1] - vin[ik]
        k[ik] = - vin[ik]/dvi
        assert np.abs((1-k[ik])*vin[ik] + k[ik]*vin[ik+1]) < 1e-5
        
    return k

def zero_hit_vec(vin,test_monotonicity=False,test_sc=True):
    assert vin.ndim == 1
    
    dv = vin[1:] - vin[:-1]
    if test_monotonicity: assert np.all(dv>0) or np.all(dv<0)
    k = - vin[:-1] / dv
    
    if test_sc:
        sum_hit = np.sum( (k >= 0)*(k<= 1) )
        assert sum_hit == 1 or sum_hit == 0
    
    return k
    
def zero_hit_mat(Vin,trim=True,test_monotonicity=False,test_sc=True,return_loc=False):
    # for a vector Vin it retunrs the "location" k of zero between coordiates
    # i and i+1. Precisely, it returns vector k where k[i] is such that 
    # x[i]*(1-k[i]) + x[i+1]*k[i] = 0. 
    # Size of k is size of Vin-1. 
    # If Vin is array, it performs the same thing operation over the LAST 
    # dimension of matrix (i.e. for matrix -- it works over column dimension)
    # 
    # If trim = True:
    #    k[i] are trimed to be between 0 and 1, so if k[i] is 1 this indicates
    #    that both x[i] and x[i+1] are (weakly) negative, 
    #    if k[i] is 0 then both are positive, so k[i] is between 0 and 1 only
    #    in case when x[i] and x[i+1] are of the opposite sign
    # It also can test monotonicity and single crossing (uniquness of hitting 
    # zero). Monotonicity guarantees single crossing.
    # F
    
    dV = np.diff(Vin,axis=-1)    
    if test_monotonicity: assert np.all(dV>0)
    
    k = - Vin[...,:-1] / dV
    
    if test_sc:
        sum_hit = np.sum( (k>=0)*(k<=1), axis=-1)
        assert np.all( (sum_hit==1) | (sum_hit==0) )
        
    if trim: k = np.maximum(np.minimum(k,1.0),0.0)
        
    if not return_loc:
        return k
    else:
        
        # now we have to figure out the location of k that is between 0 and 1
        k_01 = (k>=0)*(k<=1)
        # this requires single-crossing
        # but if violated this returns the location of the first crossing
        all_neg = np.all(k==1)
        all_pos = np.all(k==0)
        loc = np.argmax(k_01,axis=-1)
        loc[all_neg] = k.shape[-1]+1 # note that this is 2 more than largest        
        loc[all_pos] = -1 
        
        return k, loc


if __name__ == '__main__':
    vin = np.array([-2,-0.5, -0.1, 0.2, 0.5 , 2, 4])
    
    vin2 = np.array([4,3,2.5,1,0.5,-0.2,-0.5])
    
    v = np.array([vin,vin2]).reshape(2,1,7)
    
    assert np.all( np.abs(zero_hit_vec(vin2) - zero_hit_loop(vin2))<1e-5)
    print(zero_hit_mat(v))