#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
#import dill as pickle

#from ren_mar_pareto import v_mar_igrid, v_no_mar
#from ren_mar_alt import v_mar_igrid, v_no_marù
from marriage import v_mar_igrid ,v_no_mar
    



def ev_single(setup,V,sown,female,t,trim_lvl=0.001,decc=None):
    
    # # expected value of single person meeting a partner with a chance pmeet
    # if (female==False) & (t<2):
    #     pmeet = 0  #For consistency with the age gap between spouses men show up in mm 2 years later
    # else:
    pmeet = setup.dtype( setup.pars['pmeet_t'][t] )
    
    skip_mar = (pmeet < 1e-5)
    
    
    # do test here
    #ev_single_meet_test(setup,V,sown,female,t,
     #                                 skip_mar=skip_mar,trim_lvl=trim_lvl)
    
    EV_meet, dec = ev_single_meet(setup,V,sown,female,t,
                                      skip_mar=skip_mar,trim_lvl=trim_lvl,dec_c=decc)
    
    
    
    if female:
        M = setup.exogrid.zf_t_mat[t].T
        EV_nomeet =  np.dot(V['Female, single']['V'],M).astype(setup.dtype)
    else:
        M = setup.exogrid.zm_t_mat[t].T
        EV_nomeet =  np.dot(V['Male, single']['V'],M).astype(setup.dtype)
    
    assert EV_nomeet.dtype == setup.dtype
    assert EV_meet.dtype   == setup.dtype
    
    
    return (1-pmeet)*EV_nomeet + pmeet*EV_meet, dec
    

def ev_single_meet(setup,V,sown,female,t,skip_mar=False,trim_lvl=0.000001,dec_c=None):
    # computes expected value of single person meeting a partner
    
    # this creates potential partners and integrates over them
    # this also removes unlikely combinations of future z and partner's 
    # characteristics so we have to do less bargaining
    
    nexo = setup.pars['nexo_t'][t]
    ns = sown.size
    
    
    p_mat = setup.part_mats['Female, single'][t].T if female else setup.part_mats['Male, single'][t].T
    p_mat = p_mat.astype(setup.dtype,copy=False)
        
    V_next = np.ones((ns,nexo),dtype=setup.dtype)*(-1e10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    
    
    EV = setup.dtype(0.0)
    
    i_assets_c, p_assets_c = setup.i_a_mat[female], setup.prob_a_mat[female]
    
    npart = i_assets_c.shape[1]
    
    
    matches = setup.matches['Female, single'][t] if female else setup.matches['Male, single'][t]
    
    
    dec = np.zeros(matches['iexo'].shape,dtype=np.bool)
    morc = np.zeros(matches['iexo'].shape,dtype=np.bool)
    tht = -1*np.ones(matches['iexo'].shape,dtype=np.int32)
    iconv = matches['iconv']
    
    temp1=np.linspace(1,setup.pars['nexo_t'][t],setup.pars['nexo_t'][t])-1
    izfa=setup.all_indices(t,temp1)[1]
    izma=setup.all_indices(t,temp1)[2]
    psia=setup.all_indices(t,temp1)[3]
    index2=np.array(setup.all_indices(t,(izma,izfa,psia))[0],dtype=np.int16)
    
   
    
    
    for i in range(npart):
        if not skip_mar:
            
            # try marriage
            res = v_mar_igrid(setup,t,V,i_assets_c[:,i],inds,
                                      female=female,marriage=True)
         
        else:
            # try marriage
            res = v_no_mar(setup,t,V,i_assets_c[:,i],inds,
                                      female=female,marriage=True)
            
        
    
      
        (vfout,vmout), npr, dect, thtt, i_mar = res['Values'], res['NBS'], res['Decision'], res['theta'], res['i_mar']
        
       
        if female:
            vout = vfout
        else:
            vout = vmout
            
        dec[:,:,iconv[:,i]] = dect[:,None,:]
        tht[:,:,iconv[:,i]] = thtt[:,None,:]
        morc[:,:,iconv[:,i]] = i_mar[:,None,:]
            
    
            
        assert vout.dtype == setup.dtype
        
 

            
        V_next[:,inds] = vout
        
        EV += (p_assets_c[:,i][:,None])*np.dot(V_next,p_mat)
    
    assert EV.dtype == setup.dtype
    
    mout = matches.copy()
    mout['Decision'] = dec
    mout['M or C'] = morc
    mout['theta'] = tht
    
   
    
    return EV, mout


