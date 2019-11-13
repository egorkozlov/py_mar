#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
import dill as pickle
#from mc_tools import int_prob

from ren_mar import v_mar,v_mar2
    



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
    sig_a_partner = setup.pars['sig_partner_a']
    ns = sown.size
    eps_a_partner = setup.integration['nodes_couple'][:,0]
    npart = setup.integration['num_partners']    
    
    
    p_mat = setup.part_mats['Female, single'][t].T if female else setup.part_mats['Male, single'][t].T
   
        
    V_next = np.ones((ns,nexo))*(-1e10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    EV = 0.0
    
    for eps in eps_a_partner:
        
        spart = sown*np.exp(sig_a_partner*eps)
        if female: # TODO: this can be done better with keywords
           
            voutm,nprm= v_mar2(setup,V,True,sown,spart,inds,interpolate=False)[0:2]
            voutc,nprc = v_mar2(setup,V,False,sown,spart,inds,interpolate=False)[0:2]
            
            #Cohabitation-Marriage Choice 
            i_mar = (nprm>=nprc)
            
            
            
            
            with open('name_model.pkl', 'wb') as file:
                pickle.dump((i_mar,nprm,nprc), file)  
            vout = i_mar*voutm[0] + (1.0-i_mar)*voutc[0]
           
        else:
            
            voutm,nprm= v_mar2(setup,V,True,spart,sown,inds,interpolate=False)[0:2]
            voutc,nprc= v_mar2(setup,V,False,spart,sown,inds,interpolate=False)[0:2]
            
            #Cohabitation-Marriage Choice            
            i_mar = (nprm>=nprc)
            vout = i_mar*voutm[1] + (1.0-i_mar)*voutc[1]
            
        V_next[:,inds] = vout
        
        EV += (1/npart)*np.dot(V_next,p_mat)

    
    return EV#, p_mat, VF_next