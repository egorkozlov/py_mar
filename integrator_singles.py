#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
import dill as pickle
#from mc_tools import int_prob

from ren_mar import v_mar,v_mar2#, v_mar_igrid
from ren_mar_alt import v_mar_igrid
    



def ev_single(setup,V,sown,female,t,trim_lvl=0.001):
    # expected value of single person meeting a partner with a chance pmeet
    pmeet = setup.pars['pmeet']
    
    EV_meet = ev_single_meet2(setup,V,sown,female,t,trim_lvl=trim_lvl)
    
    if female:
        M = setup.exogrid.zf_t_mat[t].T
        EV_nomeet =  np.dot(V['Female, single']['V'],M)
    else:
        M = setup.exogrid.zm_t_mat[t].T
        EV_nomeet =  np.dot(V['Male, single']['V'],M)
    
    return (1-pmeet)*EV_nomeet + pmeet*EV_meet
    

def ev_single_meet(setup,V,sown,female,t,trim_lvl=0.001):
    # computes expected value of single person meeting a partner
    
    # this creates potential partners and integrates over them
    # this also removes unlikely combinations of future z and partner's 
    # characteristics so we have to do less bargaining
    
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
    
    for eps in eps_a_partner:
        
        spart = sown*np.exp(sig_a_partner*eps)
        if female: # TODO: this can be done better with keywords
            #TODO vmar_choice here
            vout = v_mar(setup,V,sown,spart,inds,interpolate=False)[i_vnext]
            
        else:
            vout = v_mar(setup,V,spart,sown,inds,interpolate=False)[i_vnext]
            
        V_next[:,inds] = vout
        
        EV += (1/npart)*np.dot(V_next,p_mat)
        
    
    return EV#, p_mat, VF_next

def ev_single_meet2(setup,V,sown,female,t,trim_lvl=0.001):
    # computes expected value of single person meeting a partner
    
    # this creates potential partners and integrates over them
    # this also removes unlikely combinations of future z and partner's 
    # characteristics so we have to do less bargaining
    
    nexo = setup.pars['nexo']
    ns = sown.size
    
    p_mat = setup.part_mats['Female, single'][t].T if female else setup.part_mats['Male, single'][t].T
   
        
    V_next = np.ones((ns,nexo))*(-1e10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
   
    
    EV = 0.0
    
    i_assets_c, p_assets_c = setup.i_a_mat, setup.prob_a_mat
    
    npart = i_assets_c.shape[1]
    
    
    for i in range(npart):
        
        res_m = v_mar_igrid(setup,V,i_assets_c[:,i],inds,
                                 female=female,marriage=True)
        
        (vfoutm,vmoutm), nprm = res_m['Values'], res_m['NBS']
        
        res_c = v_mar_igrid(setup,V,i_assets_c[:,i],inds,
                                 female=female,marriage=True)
        
        
        (vfoutc,vmoutc), nprc =  res_c['Values'], res_c['NBS']
        
        i_mar = (nprm>=nprc) 
        
        if female:
            vout = i_mar*vfoutm + (1-i_mar)*vfoutc
        else:
            vout = i_mar*vmoutm + (1-i_mar)*vmoutc
            
        V_next[:,inds] = vout
        
        EV += (p_assets_c[:,i][:,None])*np.dot(V_next,p_mat)

    
    return EV#, p_mat, VF_next