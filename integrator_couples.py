#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is integrator for couples

"""

import numpy as np
#from renegotiation import v_last_period_renegotiated, v_renegotiated_loop
#from ren_mar_pareto import v_ren_new, v_no_ren
#from ren_mar_alt import v_ren_new, v_no_ren
from renegotiation_unilateral import v_no_ren   
from renegotiation_unilateral import v_ren_uni
from renegotiation_bilateral import v_ren_bil
from renegotiation_vartheta import v_ren_vt
from renegotiation_decisions import v_ren_vt as v_ren_dec 
#from ren_mar_pareto import v_ren_new as ren_pareto

def ev_couple_m_c(setup,Vpostren,t,marriage,use_sparse=True,draw=False):
    # computes expected value of couple entering the next period with an option
    # to renegotiate or to break up
    
    can_divorce = setup.pars['can divorce'][t]
    if can_divorce:
        uni_div = setup.div_costs.unilateral_divorce if marriage else setup.sep_costs.unilateral_divorce
        title = setup.div_costs.eq_split if marriage else setup.sep_costs.eq_split 
        if uni_div:
            
            # choose your fighter
            if title>0.5:
                out = v_ren_vt(setup,Vpostren,marriage,t)#v_ren_uni(setup,Vpostren,marriage,t)
            else:
                out = v_ren_dec(setup,Vpostren,marriage,t)            
        else:
            out = v_ren_bil(setup,Vpostren,marriage,t)
    else:
        out = v_no_ren(setup,Vpostren,marriage,t)
           
       
    _Vren2 =out['Values'] if draw else out.pop('Values') 
    #_Vren2=out['Values']
    

    
    ###########################################################
    ###########################################################
    
    
    if len(out['Decision'].shape)>3:
        
        tk = lambda x : x[:,:,setup.theta_orig_on_fine,:]
        
        Vren = {'M':{'VR':_Vren2[0],'VC':_Vren2[1], 'VF':_Vren2[2],'VM':_Vren2[3]},
        'SF':Vpostren['Female, single'],
        'SM':Vpostren['Male, single']}
        
        # accounts for exogenous transitions
        EVr, EVc, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,use_sparse,down=False)
        
      
        if title>0.5:
            EVr=EVr[:,:,:,:,int((len(setup.ashare)-1)/2)]
            EVc=EVc[:,:,:,:,int((len(setup.ashare)-1)/2)]
            EVf=EVf[:,:,:,:,int((len(setup.ashare)-1)/2)]
            EVm=EVm[:,:,:,:,int((len(setup.ashare)-1)/2)]
      

    else:
        
        tk = lambda x : x[:,:,setup.theta_orig_on_fine]
        
        Vren = {'M':{'VR':tk(_Vren2[0]),'VC':tk(_Vren2[1]), 'VF':tk(_Vren2[2]),'VM':tk(_Vren2[3])},
            'SF':Vpostren['Female, single'],
            'SM':Vpostren['Male, single']}
        EVr, EVc, EVf, EVm = ev_couple_exo1(setup,Vren['M'],t,use_sparse,down=False)
        assert EVr.dtype == setup.dtype
        dec = out.copy()
    
        return (EVr, EVc, EVf, EVm), dec

    assert EVr.dtype == setup.dtype
    dec = out.copy()
    
    tk = lambda x : x[:,:,setup.theta_orig_on_fine,:,:]
    
    return (tk(EVr), tk(EVc), tk(EVf), tk(EVm)), dec#(EVr, EVc, EVf, EVm), dec#


def ev_couple_exo(setup,Vren,t,use_sparse=True,down=False):
    
 
    # this does dot product along 3rd dimension
    # this takes V that already accounts for renegotiation (so that is e
    # expected pre-negotiation V) and takes expectations wrt exogenous shocks
    
    # floating point math is quirky and can change dtypes occasionally
    def mmult(a,b):
        if use_sparse:
            return (a*b).astype(a.dtype,copy=False)
        else:
            return np.dot(a,b.T).astype(a.dtype,copy=False)
        
    
    nl = len(setup.exogrid.all_t_mat_by_l_spt)
    
    na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta_fine
    
    
    Vr, Vc, Vf, Vm = Vren['VR'], Vren['VC'], Vren['VF'], Vren['VM']
    EVr, EVc, EVf, EVm = np.zeros((4,na,nexo,ntheta,nl,len(setup.ashare)),dtype=setup.dtype)
    
    
    for il in range(nl):
        
        M = setup.exogrid.all_t_mat_by_l_spt[il][t] if use_sparse else setup.exogrid.all_t_mat_by_l[il][t]
        
        
        
        for itheta in range(ntheta):
            for j in range(len(setup.ashare)):
                #for ia in range(EVr.shape[0]):
                EVr[...,itheta,il,j]  = mmult(Vr[...,itheta,j],M)
                EVc[...,itheta,il,j]  = mmult(Vc[...,itheta,j],M)
                EVf[...,itheta,il,j]  = mmult(Vf[...,itheta,j],M)         
                EVm[...,itheta,il,j]  = mmult(Vm[...,itheta,j],M)            
                

    #assert not np.allclose( EV[...,0], EV[...,1])
    
    
    return EVr, EVc, EVf, EVm


def ev_couple_exo1(setup,Vren,t,use_sparse=True,down=False):
    
 
    # this does dot product along 3rd dimension
    # this takes V that already accounts for renegotiation (so that is e
    # expected pre-negotiation V) and takes expectations wrt exogenous shocks
    
    # floating point math is quirky and can change dtypes occasionally
    def mmult(a,b):
        if use_sparse:
            return (a*b).astype(a.dtype,copy=False)
        else:
            return np.dot(a,b.T).astype(a.dtype,copy=False)
        
    
    nl = len(setup.exogrid.all_t_mat_by_l_spt)
    
    na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta
    
    
    Vr, Vc, Vf, Vm = Vren['VR'], Vren['VC'], Vren['VF'], Vren['VM']
    EVr, EVc, EVf, EVm = np.zeros((4,na,nexo,ntheta,nl),dtype=setup.dtype)
    
    
    for il in range(nl):
        
        M = setup.exogrid.all_t_mat_by_l_spt[il][t] if use_sparse else setup.exogrid.all_t_mat_by_l[il][t]
        
        
        
        for itheta in range(ntheta):
            EVr[...,itheta,il]  = mmult(Vr[...,itheta],M)
            EVc[...,itheta,il]  = mmult(Vc[...,itheta],M)
            EVf[...,itheta,il]  = mmult(Vf[...,itheta],M)         
            EVm[...,itheta,il]  = mmult(Vm[...,itheta],M)            
            

    #assert not np.allclose( EV[...,0], EV[...,1])
    
    
    return EVr, EVc, EVf, EVm




