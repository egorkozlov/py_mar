#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is integrator for couples

"""

import numpy as np
#from renegotiation import v_last_period_renegotiated, v_renegotiated_loop
from ren_mar_alt import v_ren_new
    

def ev_couple_m_c(setup,Vpostren,t,marriage,use_sparse=True):
    # computes expected value of couple entering the next period with an option
    # to renegotiate or to break up
    
    
   
    out = v_ren_new(setup,Vpostren,marriage,t)
    #_Vren2 = out.pop('Values') # THIS REMOVES VALUES FROM OUT ASSUMING WE DO NOT NEED THEM
    _Vren2=out['Values']
    dec = out
    
    
    tk = lambda x : x[:,:,setup.theta_orig_on_fine]
    
    Vren = {'M':{'V':tk(_Vren2[0]),'VF':tk(_Vren2[1]),'VM':tk(_Vren2[2])},
            'SF':Vpostren['Female, single'],
            'SM':Vpostren['Male, single']}

    
    # accounts for exogenous transitions
    
    EV, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,use_sparse)
    
    
    return (EV, EVf, EVm), dec


def ev_couple_exo(setup,Vren,t,use_sparse=True):
    
 
    # this does dot product along 3rd dimension
    # this takes V that already accounts for renegotiation (so that is e
    # expected pre-negotiation V) and takes expectations wrt exogenous shocks
    
    
    if use_sparse:
        
        M = setup.exogrid.all_t_mat_sparse_T[t]      
        def integrate_array(x):
            xout = np.zeros_like(x)
            for itheta in range(x.shape[2]):
                xout[...,itheta] = x[...,itheta]*M                
            return xout
    else:
        
        M = setup.exogrid.all_t_mat[t].T        
        def integrate_array(x):
            xout = np.zeros_like(x)
            for itheta in range(x.shape[2]):
                xout[...,itheta] = np.dot(x[...,itheta], M)                
            return xout  
        
        
        
    
    EV, EVf, EVm = tuple( (integrate_array(a) for a in (Vren['V'],Vren['VF'],Vren['VM'])) )
    
    
    return EV, EVf, EVm