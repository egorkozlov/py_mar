#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is solver for those who are couples at period 0
"""
import numpy as np
from timeit import default_timer


from optimizers import v_optimize_couple

from platform import system

if system() != 'Darwin' and system() != 'Windows':    
    nbatch_def = 500    
else:    
    nbatch_def = 17

def v_iter_couple(setup,t,EV_tuple,ushift,dec,desc,draw=True,nbatch=nbatch_def,verbose=False,
                              force_f32 = False):
    
    if verbose: start = default_timer()
    
    agrid = setup.agrid_c
    sgrid = setup.sgrid_c
    
    dtype = setup.dtype
    
    ls = setup.ls_levels
    nls = len(ls)
    
    
    # type conversion is here
    
    zf  = setup.exogrid.all_t[t][:,0]
    zm  = setup.exogrid.all_t[t][:,1]
    zftrend = setup.pars['f_wage_trend'][t]
    zmtrend = setup.pars['m_wage_trend'][t]

    psi = setup.exogrid.all_t[t][:,2]
    beta = setup.pars['beta_t'][t]
    sigma = setup.pars['crra_power']
    R = setup.pars['R_t'][t]


    
    nexo = setup.pars['nexo_t'][t]
    shp = (setup.na,nexo,setup.ntheta,len(setup.ashare))
    
    
    wf = np.exp(zf + zftrend)
    wm = np.exp(zm + zmtrend)
    
    
    dtype_here = np.float64 if force_f32 else dtype
    
    if EV_tuple is None:
        EVr_by_l, EVc_by_l, EV_fem_by_l, EV_mal_by_l = np.zeros(((4,) + shp + (nls,)), dtype=dtype )
    else:
        EVr_by_l, EVc_by_l, EV_fem_by_l, EV_mal_by_l = EV_tuple
    
    
    
    
    # type conversion
    sgrid,sigma,beta = (dtype(x) for x in (sgrid,sigma,beta))
    
    V_couple, c_opt, s_opt, x_opt = np.empty((4,)+shp,dtype)
    i_opt, il_opt = np.empty(shp,np.int16), np.empty(shp,np.int32)
    
    V_all_l = np.empty((setup.na,nexo,setup.ntheta,nls,len(setup.ashare)),dtype=dtype)
    
    theta_val = dtype(setup.thetagrid)
    
    # the original problem is max{umult*u(c) + beta*EV}
    # we need to rescale the problem to max{u(c) + beta*EV_resc}
    

    
    #Time husband contribute to build Q
    mt=1.0-setup.mlevel
    
    # this natually splits everything onto slices
    
    
    
    for j in range(len(setup.ashare)):
        #money_i = money[:,:]
      
            
        money_t = (R*agrid, wf, wm)
        EV_t = (setup.vsgrid_c,EVr_by_l[:,:,:,:,j])
        
        
        V_pure_i, c_opt_i, x_opt_i, s_opt_i, i_opt_i, il_opt_i, V_all_l_i = \
            v_optimize_couple(money_t,sgrid,EV_t,setup.mgrid,
                              setup.ucouple_precomputed_u,setup.ucouple_precomputed_x,
                              setup.ucouple_precomputed_c,
                                  ls,beta,ushift,dtype=dtype_here,mt=mt)
           
        V_ret_i = V_pure_i + psi[None,:,None]
        
        # if dtype_here != dtype type conversion happens here
        
        V_couple[:,:,:,j] = V_ret_i # this estimate of V can be improved
        c_opt[:,:,:,j] = c_opt_i 
        s_opt[:,:,:,j] = s_opt_i
        i_opt[:,:,:,j] = i_opt_i
        x_opt[:,:,:,j] = x_opt_i
        il_opt[:,:,:,j] = il_opt_i
        V_all_l[:,:,:,:,j] = V_all_l_i # we need this for l choice so it is o
    
   
    
    assert np.all(c_opt > 0)
    
    psi_r = psi[None,:,None,None].astype(setup.dtype,copy=False)
    
    # finally obtain value functions of partners
    uf, um = setup.u_part(c_opt,x_opt,il_opt,theta_val[None,None,:,None],ushift,psi_r)
    uc = setup.u_couple(c_opt,x_opt,il_opt,theta_val[None,None,:,None],ushift,psi_r)
    
    
    EVf_all, EVm_all, EVc_all  = (setup.vsgrid_c.apply_preserve_shape(x) for x in (EV_fem_by_l, EV_mal_by_l,EVc_by_l))
    

    V_all,V_fem,V_mal  =np.empty(shp,dtype),np.empty(shp,dtype),np.empty(shp,dtype)
    for j in range(len(setup.ashare)):V_all[:,:,:,j] = uc[:,:,:,j] + beta*np.take_along_axis(np.take_along_axis(EVc_all[:,:,:,:,j],i_opt[...,None,j],0),il_opt[...,None,j],3).squeeze(axis=3)
    for j in range(len(setup.ashare)):V_fem[:,:,:,j] = uf[:,:,:,j] + beta*np.take_along_axis(np.take_along_axis(EVf_all[:,:,:,:,j],i_opt[...,None,j],0),il_opt[...,None,j],3).squeeze(axis=3)
    for j in range(len(setup.ashare)):V_mal[:,:,:,j] = um[:,:,:,j] + beta*np.take_along_axis(np.take_along_axis(EVm_all[:,:,:,:,j],i_opt[...,None,j],0),il_opt[...,None,j],3).squeeze(axis=3)
   
    
      ################################
    #HERE CHOICE OF ASSETS SPLITS
    ################################
    
    sharea=V_couple.shape[-1]
    na, nexo, ntheta = setup.na, setup.pars['nexo_t'][t], setup.ntheta_fine
    Vexp=V_couple.copy()

    #Now choose the bet division of assets
    Vexp_max=np.argmax(Vexp,axis=-1)#np.argmax(Vexp,axis=-1)#
   
    dime=Vexp.shape[-1]
    
    Vexp_max_a=np.max(Vexp,axis=-1)
    #Vexp_max=np.full(Vexp_max.shape,int((len(setup.ashare)-1)/2))
    Vexp_max_b=np.reshape(np.repeat(Vexp_max_a,dime),Vexp.shape)
   
    #Correct for symmetries
    temp1=np.cumsum(np.ones(Vexp.shape,dtype=np.int32),axis=1)-1 
    index1=setup.all_indices(t,temp1)[0]
    izf=setup.all_indices(t,temp1)[1]
    izm=setup.all_indices(t,temp1)[2]
    psia=setup.all_indices(t,temp1)[3]
    index2=setup.all_indices(t,(izm,izf,psia))[0]
    Vexp2t=Vexp[:,index2[0,:,0,0],:,:]
    indt=np.linspace(setup.ntheta_fine,1,setup.ntheta_fine,dtype=np.int16)-1
    Vexptt=Vexp2t[:,:,indt,:]
    indt2=np.linspace(len(setup.ashare),1,len(setup.ashare),dtype=np.int16)-1
    Vexp2=Vexptt[:,:,:,indt2]
    
    Vexp3=Vexp[:,index2[0,:,0,0],:,:]
    
    
    Vexp_max2=np.argmax(Vexp2,axis=-1)
    Vexp_max2_a=np.max(Vexp2,axis=-1)
    Vexp_max2_b=np.reshape(np.repeat(Vexp_max2_a,dime),Vexp.shape)
    doppio=(abs(Vexp2-Vexp_max2_b)<1e-08)
    doppio2=np.array(np.cumsum(doppio,axis=-1),dtype=np.int16)[:,:,:,-1]!=1
    #Vexp_max2[doppio2]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[doppio2]
    
    Vexp_max3=np.argmax(Vexptt,axis=-1)#np.argmax(Vexp3,axis=-1)
    
    
    wha=(izf==izm)[:,:,:,0]
    whb=Vexp_max!=Vexp_max2
    
    whc=(abs(Vexp_max+Vexp_max3-len(setup.ashare)+1)>1e-12)
    
    wh=np.full(whb.shape,False,dtype=bool)# (whb) | (whc)  | (wha) #
    
    Vexp_max[wh]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[wh]
    Vexp_max2[wh]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[wh]
    Vexp_max3[wh]=np.full(Vexp_max2.shape,int((len(setup.ashare)-1)/2))[wh]
                  
 
   
    temp=np.cumsum(np.ones(Vexp.shape,dtype=np.int32),axis=-1)-1
    
    
    #Take right age using a mask
    #columns=np.linspace(1,len(Vexp[0,:]),len(Vexp[0,:]),dtype=np.int16)
    
    col=np.reshape(np.repeat(Vexp_max,dime),Vexp.shape,order='C')
    mask=(col-temp==0)

  


    #########################################################
    ##NOW RESHAPE USING AND REDUCE DIMENSION USING ARGMAX
    #####################################################
    V_all=np.reshape(V_all[mask],Vexp_max.shape)
    V_fem=np.reshape(V_fem[mask],Vexp_max.shape)
    V_mal=np.reshape(V_mal[mask],Vexp_max.shape)
    c_opt=np.reshape(c_opt[mask],Vexp_max.shape)
    s_opt=np.reshape(s_opt[mask],Vexp_max.shape)
    i_opt=np.reshape(i_opt[mask],Vexp_max.shape)
    x_opt=np.reshape(x_opt[mask],Vexp_max.shape)
    il_opt=np.reshape(il_opt[mask],Vexp_max.shape)
  
    
    v0 = np.reshape(V_all_l[:,:,:,0,:][mask],Vexp_max.shape) 
    v1 = np.reshape(V_all_l[:,:,:,1,:][mask],Vexp_max.shape) 
    V_all_l=np.concatenate((np.expand_dims(v0,axis=3),np.expand_dims(v1,axis=3)),axis=3)
    
    def r(x): return x
    
    assert V_all.dtype == dtype
    assert V_fem.dtype == dtype
    assert V_mal.dtype == dtype
    assert c_opt.dtype == dtype
    assert x_opt.dtype == dtype
    assert s_opt.dtype == dtype
  
    #assert np.allclose(V_all,np.reshape(V_couple[mask],Vexp_max.shape),atol=1e-4,rtol=1e-1)
    try:
        assert np.allclose(V_all,V_couple,atol=1e-6,rtol=1e-5)
    except:
        #print('max difference in V is {}'.format(np.max(np.abs(V_all-np.reshape(V_couple[mask],Vexp_max.shape)))))
        pass
    
    
    ################################################
    #ADJUST DECISION ACCORDING TO DIVISION OF ASSETS
    ################################################         
    dec['Decision']=np.reshape(dec['Decision'][mask],Vexp_max.shape)
    dec['thetas']=np.reshape(dec['thetas'][mask],Vexp_max.shape)
     
    if desc=='Couple, C':
        dec['Cohabitation preferred to Marriage']=np.reshape(dec['Cohabitation preferred to Marriage'][mask],Vexp_max.shape)
   
 
    dec['assdev']=Vexp_max
    print('The mean asset share is {},conditional on diovrce, conditional on div is {}'.format(np.mean(setup.ashare[Vexp_max]),np.mean(setup.ashare[Vexp_max[1:30,:,:]][dec['thetas'][1:30,:,:]==-1])))


    if draw:
        mask=np.full(mask.shape,False)
        mask[:,:,:,int((len(setup.ashare)-1)/2)]=True
 
        dec['Values']=(np.reshape(dec['Values'][0][mask],Vexp_max.shape),
                       np.reshape(dec['Values'][1][mask],Vexp_max.shape),
                       np.reshape(dec['Values'][2][mask],Vexp_max.shape),
                       np.reshape(dec['Values'][3][mask],Vexp_max.shape))
   # return r(V_all), r(V_fem), r(V_mal), r(c_opt), r(x_opt), r(s_opt), il_opt, r(V_all_l), r(EVc_all), dec
    return r(V_all), r(V_fem), r(V_mal), r(c_opt), r(x_opt), r(s_opt), il_opt, dec






