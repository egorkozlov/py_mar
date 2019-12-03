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
    
    EV, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,use_sparse,down=False)
    
    
    return (EV, EVf, EVm), dec


def ev_couple_exo(setup,Vren,t,use_sparse=True,down=False):
    
 
    # this does dot product along 3rd dimension
    # this takes V that already accounts for renegotiation (so that is e
    # expected pre-negotiation V) and takes expectations wrt exogenous shocks
    
    
    def mmult(a,b):
        if use_sparse:
            return a*b
        else:
            return np.dot(a,b.T)
        
    
    nl = len(setup.exogrid.all_t_mat_by_l_spt)
    
    
    na, nexo, ntheta = Vren['V'].shape
    
    
    V, Vf, Vm = Vren['V'], Vren['VF'], Vren['VM']
    EV, EVf, EVm = np.zeros((na,nexo,ntheta,nl)), np.zeros((na,nexo,ntheta,nl)), np.zeros((na,nexo,ntheta,nl))
    
    
    for il in range(nl):
        
        M = setup.exogrid.all_t_mat_by_l_spt[il][t] if use_sparse else setup.exogrid.all_t_mat_by_l[il][t]
        for itheta in range(ntheta):
            EV[...,itheta,il]  = mmult( V[...,itheta],M)
            EVf[...,itheta,il] = mmult(Vf[...,itheta],M)             
            EVm[...,itheta,il] = mmult(Vm[...,itheta],M)             
            
    
    
    
    
    return EV, EVf, EVm