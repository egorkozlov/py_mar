#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:49:58 2020

@author: egorkozlov
"""
import numpy as np
from aux_routines import first_true, last_true
from numba import njit, vectorize
from gridvec import VecOnGrid


    
def ren_bilateral_wrap(setup,vy,vfy,vmy,vfn,vmn,vf_all_s,vm_all_s,aleft_c,                       
                       ia_div_fem,ia_div_mal,rescale=True):
    # v_ren_core_interp(setup,vy,vfy,vmy,vf_n,vm_n,unilateral,show_sc=False,rescale=False)
    tgrid = setup.thetagrid_fine
    
    
    vout, vfout, vmout, thetaout, yes, ithetaout, bribe, iaout_f, iaout_m = \
        ren_loop_bilateral(vy,vfy,vmy,vfn,vmn,
                           vf_all_s,vm_all_s,aleft_c,
                           ia_div_fem,ia_div_mal,
                           setup.agrid_s,
                           tgrid)
    
    if np.any(bribe):
        print('Bribing happens in {}% of divorces'.format(round(100*np.mean(~yes & bribe)/np.mean(~yes))))
    
    def r(x): return x.astype(np.float32)        
    
    
    return {'Decision': yes, 'thetas': ithetaout,
            'Values': (r(vout), r(vfout), r(vmout)),'Divorce':(vfn,vmn)}
    
    
                    

@njit
def ren_loop_bilateral(vy,vfy,vmy,vfn,vmn,vfn_as,vmn_as,aleft_c,ia_f_def_s,ia_m_def_s,agrid_s,thtgrid):
    #print('bilateral hi!')


    #vfn = vfn
    #vmn = vmn
    
    sf = vfy - vfn.reshape(vfn.shape+(1,))
    sm = vmy - vmn.reshape(vmn.shape+(1,))
    
    na, nexo, nt = vy.shape
    
    vout = vy.copy()
    vfout = vfy.copy()
    vmout = vmy.copy()
    
    yes = np.zeros((na,nexo,nt),dtype=np.bool_)
    bribe = np.zeros((na,nexo,nt),dtype=np.bool_)
    
    
    thetaout = -1*np.ones(vout.shape,dtype=np.float32)
    ithetaout = -1*np.ones(vout.shape,dtype=np.int16)
    
    iaout_f = -1*np.ones(vout.shape,dtype=np.int16)
    iaout_m = -1*np.ones(vout.shape,dtype=np.int16)
   
    
    na_s = agrid_s.size
    
    for ia in range(na):
        for ie in range(nexo):
            for it in range(nt):
                
                sf_i = sf[ia,ie,it]
                sm_i = sm[ia,ie,it]
                
                tht = thtgrid[it]
                
                
                if sf_i >= 0 and sm_i >= 0:
                    yes[ia,ie,it] = True
                    thetaout[ia,ie,it] = tht
                    ithetaout[ia,ie,it] = it
                    continue
                else:
                    
                    vout_div_def = tht*vfn[ia,ie] +  (1-tht)*vmn[ia,ie]
                    vfout_div_def = vfn[ia,ie]
                    vmout_div_def = vmn[ia,ie]
                    
                    if sf_i < 0 and sm_i < 0:                        
                        vout[ia,ie,it] = vout_div_def
                        vfout[ia,ie,it] = vfout_div_def
                        vmout[ia,ie,it] = vmout_div_def
                        continue
                    else:
                        # if only one person wants do divorce -- possible
                        # to find reallocation of assets such that both could
                        # agree.
                        ia_m_def = ia_m_def_s[ia,ie,it]
                        ia_f_def = ia_f_def_s[ia,ie,it]
                        
                        ia_m_new = ia_m_def
                        ia_f_new = ia_f_def
                        
                        a_left = aleft_c[ia,ie,it]
                        
                        
                        do_divorce = False
                        found = False
                        

                        if sf_i < 0 and sm_i > 0:
                            # f bribes out
                            #print('at point {} m bribes out'.format((ia,ie,it)))
                            # increase ia_m
                            for ia_m_new in range(ia_m_def+1,na_s):
                                if agrid_s[ia_m_new] > a_left:
                                    break
                                
                                found = False
                                for ia_f_new in range(ia_f_new,-1,-1):
                                    if agrid_s[ia_f_new] + agrid_s[ia_m_new] <= a_left:
                                        found=True
                                        break
                                    
                                if found:
                                    sf_i_new  = vfy[ia,ie,it] - vfn_as[ia_f_new,ie]
                                    sm_i_new  = vmy[ia,ie,it] - vmn_as[ia_m_new,ie]
                                    if sf_i_new < 0 and sm_i_new < 0:
                                        do_divorce = True
                                        #print('divorce happens: f bribes out!')
                                        #print((ia_m_def,sm_i,ia_f_def,sf_i))
                                        #print((ia_f_new,sf_i_new,ia_m_new,sm_i_new))                                       
                                        break
                        
                        
                        if sm_i < 0 and sf_i > 0:
                            # m bribes out
                            # increase ia_f
                            for ia_f_new in range(ia_f_def+1,na_s):
                                if agrid_s[ia_f_new] > a_left:
                                    break
                                
                                found = False
                                for ia_m_new in range(ia_m_new,-1,-1):
                                    if agrid_s[ia_m_new] + agrid_s[ia_f_new] <= a_left:
                                        found=True
                                        break
                                    
                                if found:
                                    sf_i_new  = vfy[ia,ie,it] - vfn_as[ia_f_new,ie]
                                    sm_i_new  = vmy[ia,ie,it] - vmn_as[ia_m_new,ie]
                                    if sf_i_new < 0 and sm_i_new < 0:
                                        do_divorce = True
                                        #print('divorce happens: m bribes out!')
                                        #print((ia_m_def,sm_i,ia_f_def,sf_i))
                                        #print((ia_f_new,sf_i_new,ia_m_new,sm_i_new))
                                        break
                        
                                
                        if not do_divorce:
                            yes[ia,ie,it] = True
                            thetaout[ia,ie,it] = thtgrid[it]
                            ithetaout[ia,ie,it] = it
                            continue
                        
                        # else we do_divorce   
                        assert found
                        bribe[ia,ie,it] = True
                        vfout[ia,ie,it] = vfn_as[ia_f_new,ie]
                        vmout[ia,ie,it] = vmn_as[ia_m_new,ie]
                        vout[ia,ie,it] = tht*vfn_as[ia_f_new,ie] + \
                                        (1-tht)*vmn_as[ia_m_new,ie]
                        
                        iaout_f[ia,ie,it] = ia_f_new
                        iaout_m[ia,ie,it] = ia_m_new
                        
                        
                        continue
                        
            
    return vout, vfout, vmout, thetaout, yes, ithetaout, bribe, iaout_f, iaout_m
    