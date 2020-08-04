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

    



    ############################################################
    #Now, if title based regime get the optimal division rule
    ###########################################################
    if len(out['Values'][0].shape)>3:
        
        def mmult(a,b):
            if use_sparse:
                return (a*b).astype(a.dtype,copy=False)
            else:
                return np.dot(a,b.T).astype(a.dtype,copy=False)
            
            
        #Get expected value after prod shocks
        M = setup.exogrid.all_t_mat_psi_spt[t] if use_sparse else setup.exogrid.all_t_mat_psi[t]
            
        sharea=out['Values'][0].shape[-1]
        na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta_fine
        Vexp = np.zeros((na,nexo,ntheta,sharea),dtype=setup.dtype)
        for itheta in range(setup.ntheta_fine):
            for ishare in range(sharea):
                Vexp[...,itheta,ishare]=mmult(out['Values'][0][...,itheta,ishare],M)
                
        #Now choose the bet division of assets
        Vexp_max=np.argmax(Vexp,axis=-1)#np.argmax(out['Values'][1],axis=-1)#
        temp=np.cumsum(np.ones(out['Values'][0].shape,dtype=np.int32),axis=-1)-1
        
        #Take right age using a mask
        #columns=np.linspace(1,len(Vexp[0,:]),len(Vexp[0,:]),dtype=np.int16)
        dime=out['Values'][1].shape[-1]
        col=np.reshape(np.repeat(Vexp_max,dime),out['Values'][0].shape)
        mask=(col-temp==0)
    
        out['Values']=(np.reshape(out['Values'][0][mask],Vexp_max.shape),
                       np.reshape(out['Values'][1][mask],Vexp_max.shape),
                       np.reshape(out['Values'][2][mask],Vexp_max.shape),
                       np.reshape(out['Values'][3][mask],Vexp_max.shape))
        
        out['Decision']=np.reshape(out['Decision'][mask],Vexp_max.shape)
        out['thetas']=np.reshape(out['thetas'][mask],Vexp_max.shape)
        out['Divorce']=(np.reshape(out['Divorce'][0][mask[:,:,10,:]],(Vexp_max.shape[0],Vexp_max.shape[1])),
                        np.reshape(out['Divorce'][1][mask[:,:,10,:]],(Vexp_max.shape[0],Vexp_max.shape[1])))
        
        
        if not marriage:
            out['Cohabitation preferred to Marriage']=np.reshape(out['Cohabitation preferred to Marriage'][mask],Vexp_max.shape)
       
        #Get the preferred division of assets
        out['assdev']=Vexp_max
        print('The mean asset share is {},conditional on diovrce, conditional on div is {}'.format(np.mean(setup.ashare[Vexp_max]),np.mean(setup.ashare[Vexp_max][out['thetas']==-1])))
        
        
       
    _Vren2 =out['Values'] if draw else out.pop('Values') 
    #_Vren2=out['Values']
    dec = out
    
    
    tk = lambda x : x[:,:,setup.theta_orig_on_fine]
    
    Vren = {'M':{'VR':tk(_Vren2[0]),'VC':tk(_Vren2[1]), 'VF':tk(_Vren2[2]),'VM':tk(_Vren2[3])},
            'SF':Vpostren['Female, single'],
            'SM':Vpostren['Male, single']}
    
    ###########################################################
    ###########################################################
    
    # accounts for exogenous transitions
    EVr, EVc, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,use_sparse,down=False)
    
    
    assert EVr.dtype == setup.dtype
    
    return (EVr, EVc, EVf, EVm), dec


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