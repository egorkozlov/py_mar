#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is integrator for couples

"""

import numpy as np

def ev_couple_after_savings(setup,Vren,use_sparse=True):
    # takes renegotiated values as input
    
    # this does dot product along 3rd dimension
    
    
    
    
    
    
    if use_sparse:
        
        M = setup.exogrid.all_t_mat_sparse_T[0]      
        def integrate_array(x):
            xout = np.zeros_like(x)
            for itheta in range(x.shape[2]):
                xout[...,itheta] = x[...,itheta]*M                
            return xout
    else:
        
        M = setup.exogrid.all_t_mat[0].T        
        def integrate_array(x):
            xout = np.zeros_like(x)
            for itheta in range(x.shape[2]):
                xout[...,itheta] = np.dot(x[...,itheta], M)                
            return xout  
        
        
        
    
    EV, EVf, EVm = tuple( (integrate_array(a) for a in (Vren['V'],Vren['VF'],Vren['VM'])) )
    
    
    return EV, EVf, EVm