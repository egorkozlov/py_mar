#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
#from mc_tools import int_prob

from ren_mar import v_mar
    

def ev_after_savings_grid_all_z(setup,V,sown,female,t,trim_lvl=0.01):
    # this takes gender as argument so should be called twice
    
    nexo = setup.pars['nexo']
    sig_a_partner = setup.pars['sig_partner_a']
    ns = sown.size
    eps_a_partner = setup.integration['nodes_couple'][:,0]
    npart = setup.integration['num_partners']    
    
    
    p_mat = setup.part_mats['Female, single'][t].T if female else setup.part_mats['Male, single'][t].T
    i_vnext = 0 if female else 1
        
    V_next = np.ones((ns,nexo))*(-1e10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    EV = 0.0
    
    for ipart in range(npart):
        sm = sown*np.exp(sig_a_partner*eps_a_partner[ipart])
        #vout2 = v_after_mar_grid(setup,V,sown,sm,inds)[i_vnext]        
        vout = v_mar(setup,V,sown,sm,inds,interpolate=False)[i_vnext]
        V_next[:,inds] = vout
        #print(np.max(np.abs(vout-vout2)))
        
        EV += (1/npart)*np.dot(V_next,p_mat)
        
    
    return EV#, p_mat, VF_next