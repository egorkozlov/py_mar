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
                out = v_ren_dec(setup,Vpostren,marriage,t)#v_ren_uni(setup,Vpostren,marriage,t)
            else:
                out = v_ren_dec(setup,Vpostren,marriage,t)            
        else:
            out = v_ren_bil(setup,Vpostren,marriage,t)
    else:
        out = v_no_ren(setup,Vpostren,marriage,t)
           
       
    _Vren2 =out['Values'] if draw else out.pop('Values') 
    #_Vren2=out['Values']
    
    
    
    tk = lambda x : x[:,:,setup.theta_orig_on_fine]
    
    Vren = {'M':{'VR':tk(_Vren2[0]),'VC':tk(_Vren2[1]), 'VF':tk(_Vren2[2]),'VM':tk(_Vren2[3])},
            'SF':Vpostren['Female, single'],
            'SM':Vpostren['Male, single']}
    
    ###########################################################
    ###########################################################
    
    
    if len(out['Values'][1].shape)>3:
        # accounts for exogenous transitions
        EVr, EVc, EVf, EVm = ev_couple_exo(setup,Vren['M'],t,use_sparse,down=False)
        
        
        
        sharea=EVr.shape[-1]
        na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta_fine
        Vexp=EVr
    
        #Now choose the bet division of assets
        Vexp_max=np.argmax(Vexp,axis=-1)#np.argmax(Vexp,axis=-1)#
       
        dime=Vexp.shape[-1]
        
        Vexp_max_a=np.max(Vexp,axis=-1)
        Vexp_max_b=np.reshape(np.repeat(Vexp_max_a,dime),Vexp.shape)
       
        #Correct for symmetries
        temp1=np.cumsum(np.ones(Vexp.shape,dtype=np.int32),axis=1)-1 
        index1=setup.all_indices(t,temp1)[0]
        izf=setup.all_indices(t,temp1)[1]
        izm=setup.all_indices(t,temp1)[2]
        psi=setup.all_indices(t,temp1)[3]
        index2=setup.all_indices(t,(izm,izf,psi))[0]
        Vexp2t=Vexp[:,index2[0,:,0,0,0],:,:,:]
        indt=np.linspace(setup.ntheta_fine,1,setup.ntheta_fine,dtype=np.int16)-1
        Vexptt=Vexp2t[:,:,indt,:,:]
        indt2=np.linspace(len(setup.ashare),1,len(setup.ashare),dtype=np.int16)-1
        Vexp2=Vexptt[:,:,:,:,indt2]
        
        Vexp3=Vexp[:,index2[0,:,0,0,0],:,:,:]
        
        
        Vexp_max2=np.argmax(Vexp2,axis=-1)
        Vexp_max2_a=np.max(Vexp2,axis=-1)
        Vexp_max2_b=np.reshape(np.repeat(Vexp_max2_a,dime),Vexp.shape)
        doppio=(abs(Vexp2-Vexp_max2_b)<1e-08)
        doppio2=np.array(np.cumsum(doppio,axis=-1),dtype=np.int16)[:,:,:,:,-1]!=1
        #Vexp_max2[doppio2]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[doppio2]
        
        Vexp_max3=np.argmax(Vexptt,axis=-1)#np.argmax(Vexp3,axis=-1)
        
        
        wha=(izf==izm)[:,:,:,:,0]
        whb=Vexp_max!=Vexp_max2
        
        whc=(abs(Vexp_max+Vexp_max3-len(setup.ashare)+1)>1e-12)
        
        wh= (whb) | (whc)  | (wha) #np.full(whb.shape,False,dtype=bool)#
        
        Vexp_max[wh]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[wh]
        Vexp_max2[wh]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[wh]
        Vexp_max3[wh]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[wh]
                      
        
    
        
       
        temp=np.cumsum(np.ones(Vexp.shape,dtype=np.int32),axis=-1)-1
        
        #Get the preferred division of assets
        if title>0.5:Vexp_max=np.full(Vexp_max.shape,int((len(setup.ashare)-1)/2))
           
        
        #Take right age using a mask
        #columns=np.linspace(1,len(Vexp[0,:]),len(Vexp[0,:]),dtype=np.int16)
        
        col=np.reshape(np.repeat(Vexp_max,dime),Vexp.shape,order='C')
        mask=(col-temp==0)
        
                
        if t==44:
            
            print('here')
        #assert np.all(np.max(Vexp,axis=-1)==np.reshape(Vexp[mask],Vexp_max.shape))
        out['Values']=(np.reshape(np.repeat(out['Values'][0][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0],
                       np.reshape(np.repeat(out['Values'][1][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0],
                       np.reshape(np.repeat(out['Values'][2][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0],
                       np.reshape(np.repeat(out['Values'][3][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0])
        
        out['Decision']=np.reshape(np.repeat(out['Decision'][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0]
        out['thetas']=np.reshape(np.repeat(out['thetas'][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0]
      
        EVr=np.reshape(EVr[mask],Vexp_max.shape)
        EVc=np.reshape(EVc[mask],Vexp_max.shape)
        EVf=np.reshape(EVf[mask],Vexp_max.shape)
        EVm=np.reshape(EVm[mask],Vexp_max.shape)

        # out['Divorce']=(np.reshape(out['Divorce'][0][mask[:,:,0,:]],(Vexp_max.shape[0],Vexp_max.shape[1])),
        #                 np.reshape(out['Divorce'][1][mask[:,:,0,:]],(Vexp_max.shape[0],Vexp_max.shape[1])))
        
        
        if not marriage:
            out['Cohabitation preferred to Marriage']=np.reshape(np.repeat(out['Cohabitation preferred to Marriage'][:,:,:,np.newaxis,:], 2, axis=3)[mask],Vexp_max.shape)[:,:,:,0]
      
    
        out['assdev']=Vexp_max[:,:,:,0]
        print('The mean asset share is {},conditional on diovrce, conditional on div is {}'.format(np.mean(setup.ashare[Vexp_max]),np.mean(setup.ashare[Vexp_max[1:30,:,:]][out['thetas'][1:30,:,:]==-1])))

        if t<=46:
            
            Vm=out['Values'][2]#Vpostren['Couple, M']['VM']#out['Divorce'][0]##
            Vexp=out['Values'][3]#Vpostren['Couple, M']['VF']##out['Divorce'][1]#  
            temp1=np.linspace(1,setup.pars['nexo_t'][t],setup.pars['nexo_t'][t])-1
            izf=setup.all_indices(t,temp1)[1]
            izm=setup.all_indices(t,temp1)[2]
            psi=setup.all_indices(t,temp1)[3]
            index2=np.array(setup.all_indices(t,(izm,izf,psi))[0],dtype=np.int16)
            indt=np.linspace(setup.ntheta_fine,1,setup.ntheta_fine,dtype=np.int16)-1
            Vexp2a=Vexp[:,index2,:]
            Vf=Vexp2a[:,:,indt]
            
            aaa=np.where(abs(Vm-Vf)>1e-07)
            assert np.all(abs(Vm-Vf)<1e-09)
            print('difference in integrator is {} in {}'.format(np.max(abs(Vm-Vf)),t))
            print('SINGLE DIFF is {} in {}'.format(np.max(abs(Vpostren['Male, single']['V']-Vpostren['Female, single']['V'])),t))
            
            #assert np.all(abs(Vm-Vf)<1e-07)
            #assert np.all(abs(Vpostren['Male, single']['V']-Vpostren['Female, single']['V'])<1e-07)
    
            #Eliminate below
            thetam=out['thetas']
            thetaf1=out['thetas'][:,index2,:]
            thetaf=thetaf1[:,:,indt]
           
            
            thetac=(thetam+thetaf!=setup.ntheta_fine)[thetam!=thetaf]
            
           
            print('theta mean is {} in {}'.format(np.mean(thetac),t))
            
      

    else:
        EVr, EVc, EVf, EVm = ev_couple_exo1(setup,Vren['M'],t,use_sparse,down=False)
    assert EVr.dtype == setup.dtype
    dec = out.copy()
    
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



