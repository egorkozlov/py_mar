#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
#import dill as pickle

from ren_mar_alt import v_mar_igrid, v_no_mar
    



def ev_single(setup,V,sown,female,t,trim_lvl=0.001):
    # expected value of single person meeting a partner with a chance pmeet
    pmeet = setup.pars['pmeet_t'][t]
    
    skip_mar = (pmeet < 1e-5)
    
    EV_meet, dec = ev_single_meet(setup,V,sown,female,t,
                                      skip_mar=skip_mar,trim_lvl=trim_lvl)
    
    
    
    if female:
        M = setup.exogrid.zf_t_mat[t].T
        EV_nomeet =  np.dot(V['Female, single']['V'],M)
    else:
        M = setup.exogrid.zm_t_mat[t].T
        EV_nomeet =  np.dot(V['Male, single']['V'],M)
    
    return (1-pmeet)*EV_nomeet + pmeet*EV_meet, dec
    

def ev_single_meet(setup,V,sown,female,t,skip_mar=False,trim_lvl=0.001):
    # computes expected value of single person meeting a partner
    
    # this creates potential partners and integrates over them
    # this also removes unlikely combinations of future z and partner's 
    # characteristics so we have to do less bargaining
    
    nexo = setup.pars['nexo_t'][t]
    ns = sown.size
    
    
    p_mat = setup.part_mats['Female, single'][t].T if female else setup.part_mats['Male, single'][t].T
   
        
    V_next = np.ones((ns,nexo))*(-1e10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    
    
    EV = 0.0
    
    i_assets_c, p_assets_c = setup.i_a_mat, setup.prob_a_mat
    
    npart = i_assets_c.shape[1]
    
    
    matches = setup.matches['Female, single'][t] if female else setup.matches['Male, single'][t]
    
    
    dec = np.zeros(matches['iexo'].shape,dtype=np.bool)
    morc = np.zeros(matches['iexo'].shape,dtype=np.bool)
    tht = -1*np.ones(matches['iexo'].shape,dtype=np.int32)
    iconv = matches['iconv']
    
    for i in range(npart):
        if not skip_mar:
            # try marriage
            res_m = v_mar_igrid(setup,t,V,i_assets_c[:,i],inds,
                                     female=female,marriage=True)
            
            
            res_c = v_mar_igrid(setup,t,V,i_assets_c[:,i],inds,
                                     female=female,marriage=False)
        else:
            # try marriage
            res_m = v_no_mar(setup,t,V,i_assets_c[:,i],inds,
                                     female=female,marriage=True)
            
            
            res_c = v_no_mar(setup,t,V,i_assets_c[:,i],inds,
                                     female=female,marriage=False)
        
        
        
        (vfoutm,vmoutm), nprm, decm, thtm = res_m['Values'], res_m['NBS'], res_m['Decision'], res_m['theta']
        
        # try cohabitation
        (vfoutc, vmoutc), nprc, decc, thtc =  res_c['Values'], res_c['NBS'], res_c['Decision'], res_c['theta']
        
        # choice is made based on Nash Surplus value
        i_mar = (nprm>=nprc) 
        
        if female:
            vout = i_mar*vfoutm + (1-i_mar)*vfoutc
        else:
            vout = i_mar*vmoutm + (1-i_mar)*vmoutc
            
        dec[:,:,iconv[:,i]] = (i_mar*decm + (1-i_mar)*decc)[:,None,:]
        tht[:,:,iconv[:,i]] = (i_mar*thtm + (1-i_mar)*thtc)[:,None,:]
        morc[:,:,iconv[:,i]] = i_mar[:,None,:]
            
        V_next[:,inds] = vout
        
        EV += (p_assets_c[:,i][:,None])*np.dot(V_next,p_mat)

    mout = matches.copy()
    mout['Decision'] = dec
    mout['M or C'] = morc
    mout['theta'] = tht
    
    return EV, mout