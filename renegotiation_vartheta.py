#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 19:33:42 2020

@author: egorkozlov
"""

import numpy as np
from aux_routines import first_true, last_true
from numba import njit, cuda, f4
from gridvec import VecOnGrid


def v_ren_vt(setup,V,marriage,t,return_extra=False,return_vdiv_only=False,rescale=True,
             thetafun=None):
    # this returns value functions for couple that entered the period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    # 
    # combine = True creates matrix (n_sc-by-n_inds)
    # combine = False assumed that n_sc is the same shape as n_inds and creates
    # a flat array.
     
    #Get Divorce or Separation Costs
    if marriage:
        dc = setup.div_costs
        is_unil = dc.unilateral_divorce # whether to do unilateral divorce at all
    else:
        dc = setup.sep_costs
        is_unil = dc.unilateral_divorce # whether to do unilateral divorce at all
    
    assert is_unil
    
    ind, izf, izm, ipsi = setup.all_indices(t+1)
    
    sc = setup.agrid_c
    
    #income_share_f = (np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) ).squeeze()
    
    
    if thetafun is None:
        def thetafun(tht): return tht, 1-tht
        #def thetafun(tht): return 0.5*np.ones_like(tht), 0.5*np.ones_like(tht)
            
    
    # values of divorce
    vf_n, vm_n = v_div_vartheta(
        setup, dc, t, sc,
        V['Male, single']['V'], V['Female, single']['V'],
        izf, izm, cost_fem=dc.money_lost_f, cost_mal=dc.money_lost_m, fun=thetafun)
    
    
    
    
    if return_vdiv_only:
        return {'Value of Divorce, male': vm_n,
                'Value of Divorce, female': vf_n}
    
    
    
    expnd = lambda x : setup.v_thetagrid_fine.apply(x,axis=2)
    
    if marriage:
        # if couple is married already
        v_y = expnd(V['Couple, M']['V'])
        vf_y = expnd(V['Couple, M']['VF'])
        vm_y = expnd(V['Couple, M']['VM'])
    else:
        # stay in cohabitation
        v_y_coh = expnd(V['Couple, C']['V'])
        vf_y_coh = expnd(V['Couple, C']['VF'])
        vm_y_coh = expnd(V['Couple, C']['VM'])
        # switch to marriage
        v_y_mar = expnd(V['Couple, M']['V'])
        vf_y_mar = expnd(V['Couple, M']['VF'])
        vm_y_mar = expnd(V['Couple, M']['VM'])
        # switching criterion
        #switch = (vf_y_mar>vf_y_coh) & (vm_y_mar>vm_y_coh)
        switch = (v_y_mar>= v_y_coh)
        
        v_y = switch*v_y_mar + (~switch)*v_y_coh
        vf_y = switch*vf_y_mar + (~switch)*vf_y_coh
        vm_y = switch*vm_y_mar + (~switch)*vm_y_coh
        
    vf_n, vm_n = [x.astype(v_y.dtype) for x in (vf_n,vm_n)] # type conversion
    
    v_out, vf_out, vm_out, itheta_out = \
        v_ren_core(v_y, vf_y, vm_y, vf_n, vm_n, setup.thetagrid_fine, 
                   rescale=rescale)
        
        
    v_out, vf_out, vm_out, itheta_out = \
        v_ren_gpu(v_y, vf_y, vm_y, vf_n, vm_n, setup.thetagrid_fine, 
                   rescale=rescale)
        
    def r(x): return x.astype(np.float32)
        
    result =  {'Decision': (itheta_out>=0), 'thetas': itheta_out,
                'Values': (r(v_out), r(vf_out), r(vm_out)),'Divorce':(vf_n,vm_n)}
    
    
    if not marriage:
        result['Cohabitation preferred to Marriage'] = ~switch
        
       
    extra = {'Values':result['Values'],
             'Value of Divorce, male': vm_n, 'Value of Divorce, female': vf_n}
    
    
    if not return_extra:
        return result
    else:
        return result, extra
    
    




def v_div_vartheta(setup,dc,t,sc,Vmale,Vfemale,izf,izm,
                   cost_fem=0.0,cost_mal=0.0, fun=lambda x : (x,1-x) ):
    # this produces value of divorce for gridpoints given possibly different
    # shares of how assets are divided. 
    # Returns Vf_divorce, Vm_divorce -- values of singles in case of divorce
    # matched to the gridpionts for couples
    
    # optional cost_fem and cost_mal are monetary costs of divorce
    
    from renegotiation_unilateral import v_div_allsplits
    
    shrs = setup.thetagrid_fine
    
    # these are interpolation points
    
    Vm_divorce_M, Vf_divorce_M = v_div_allsplits(setup,dc,t,sc,
                                                 Vmale,Vfemale,izm,izf,
                                shrs=shrs,cost_fem=cost_fem,cost_mal=cost_mal)
    
    # share of assets that goes to the female
    # this has many repetative values but it turns out it does not matter much
    
    
    
    
    share_fem, share_mal = fun(setup.thetagrid_fine)
    fem_gets = VecOnGrid(np.array(shrs),share_fem)
    mal_gets = VecOnGrid(np.array(shrs),share_mal)
    
    i_fem = fem_gets.i
    wn_fem = fem_gets.wnext
    
    i_mal = mal_gets.i
    wn_mal = mal_gets.wnext
    
    
    
    Vf_divorce = (1-wn_fem[None,None,:])*Vf_divorce_M[:,:,i_fem] + \
                     wn_fem[None,None,:]*Vf_divorce_M[:,:,i_fem+1]
    
    Vm_divorce = (1-wn_mal[None,None,:])*Vm_divorce_M[:,:,i_mal] + \
                     wn_mal[None,None,:]*Vm_divorce_M[:,:,i_mal+1]
                
    
                
    return Vf_divorce, Vm_divorce

@njit
def v_ren_core(v_y, vf_y, vm_y, vf_n, vm_n, thtgrid, rescale = True):
    
    na, ne, nt = v_y.shape
    
    v_out = v_y.copy()
    vm_out = vm_y.copy()
    vf_out = vf_y.copy()
    
    itheta_out = np.full(v_y.shape,-1,dtype=np.int16)
    
    i_good_current_f = (vf_y >= vf_n)
    i_good_current_m = (vm_y >= vm_n)
    
    # good -- no change
    
    i_sq = (i_good_current_f) & (i_good_current_m)
    i_anyway = (~i_good_current_f) & (~i_good_current_m)
    
    
    for ia in range(na):
        for ie in range(ne):
            for it in range(nt):
                
                vf_no = vf_n[ia,ie,it]
                vm_no = vm_n[ia,ie,it]
                
                if i_sq[ia,ie,it]:
                    # no search just fill the value
                    itheta_out[ia,ie,it] = it                    
                    continue
                    
                if i_anyway[ia,ie,it]:
                    # no search
                    tht = thtgrid[it]
                    v_out[ia,ie,it] = tht*vf_no + (1-tht)*vm_no
                    vf_out[ia,ie,it] = vf_no
                    vm_out[ia,ie,it] = vm_no
                    itheta_out[ia,ie,it] = -1
                    continue
                
                # in the points left one guy agrees and one disagrees
                
                # run two loops: forward and backward
                # see if there is anything to replace
                
                it_ren = -1
                
                found_increase = False
                found_decrease = False
                
                
                # these loops can be improved by monotonicity
                for it_increase in range(it+1,nt):
                    if (vf_y[ia,ie,it_increase] >= vf_no and vm_y[ia,ie,it_increase] >= vm_no):
                        found_increase = True
                        break
                
                for it_decrease in range(it-1,-1,-1):
                    if (vf_y[ia,ie,it_decrease] >= vf_no and vm_y[ia,ie,it_decrease] >= vm_no):
                        found_decrease = True
                        break
                    
                
                if found_increase and found_decrease:
                    dist_increase = it_increase - it
                    dist_decrease = it - it_decrease
                    
                    if dist_increase != dist_decrease:
                        it_ren = it_increase if dist_increase < dist_decrease else it_decrease
                    else:
                        # tie breaker
                        dist_mid_inc = np.abs(it_increase - (nt/2))
                        dist_mid_dec = np.abs(it_decrease - (nt/2))
                        it_ren = it_increase if dist_mid_inc < dist_mid_dec else it_decrease
                    
                elif found_increase and not found_decrease:
                    it_ren = it_increase
                elif found_decrease and not found_increase:
                    it_ren = it_decrease
                else:
                    it_ren = -1 # check this!
                    
                # finally fill the values    
                    
                if it_ren == -1:
                    tht = thtgrid[it]
                    v_out[ia,ie,it] = tht*vf_no + (1-tht)*vm_no
                    vf_out[ia,ie,it] = vf_no
                    vm_out[ia,ie,it] = vm_no
                    itheta_out[ia,ie,it] = -1
                else:
                    # here we need to rescale
                    
                    if rescale:
                        tht_old = thtgrid[it]
                        tht_new = thtgrid[it_ren]
                        factor = np.maximum( (1-tht_old)/(1-tht_new), tht_old/tht_new )
                    else:
                        factor = 1
                        
                    v_out[ia,ie,it] = factor*v_y[ia,ie,it_ren]
                    vf_out[ia,ie,it] = vf_y[ia,ie,it_ren]
                    vm_out[ia,ie,it] = vm_y[ia,ie,it_ren]
                    itheta_out[ia,ie,it] = it_ren
    
    assert np.all(vf_out >= vf_n - 1e-4)
    assert np.all(vm_out >= vm_n - 1e-4)
    
    return v_out, vf_out, vm_out, itheta_out
        

        
def v_ren_gpu(v_y, vf_y, vm_y, vf_n, vm_n, thtgrid, rescale = True):
    
                    
    na, ne, nt = v_y.shape
    assert rescale, 'no rescale is not implemented'
    
    assert v_y.shape[2] < 500 
    
    v_out = v_y.copy()
    vm_out = vm_y.copy()
    vf_out = vf_y.copy()    
    itheta_out = np.full(v_y.shape,-1,dtype=np.int16)
    
    for ia in range(na):
        for ie in range(ne):
            threadsperblock = (1, 1, nt)
                
            b_a = 1
            b_exo = 1
            b_theta = 1
            
            blockspergrid = (b_a, b_exo, b_theta)
            
            v_yi, vf_yi, vm_yi = [cuda.to_device(
                                    np.ascontiguousarray(x[ia:(ia+1),ie:(ie+1),:].copy())
                                                ) for x in (v_y, vf_y, vm_y)]
            
            vf_ni, vm_ni = [cuda.to_device(
                                            np.ascontiguousarray(x[ia:(ia+1),ie:(ie+1),:].copy())
                                          ) for x in (vf_n,vm_n)]
            
            
            v_outi, vf_outi, vm_outi, itheta_outi = [cuda.to_device(
                                     np.ascontiguousarray(x[ia:(ia+1),ie:(ie+1),:].copy())
                                    ) for x in (v_out, vm_out, vf_out, itheta_out)]
                                            
            thtgrid = cuda.to_device(thtgrid)
            
            cuda_ker[blockspergrid, threadsperblock](v_yi, vf_yi, vm_yi, vf_ni, vm_ni, 
                                            thtgrid, v_outi, vm_outi, vf_outi, itheta_outi)
        
        
            v_out[ia:(ia+1),ie:(ie+1),:] = v_outi
            vm_out[ia:(ia+1),ie:(ie+1),:] = vm_outi
            vf_out[ia:(ia+1),ie:(ie+1),:] = vf_outi
            itheta_out[ia:(ia+1),ie:(ie+1),:] = itheta_outi
            
    
    return v_out, vf_out, vm_out, itheta_out
    

@cuda.jit   
def cuda_ker(v_y, vf_y, vm_y, vf_n, vm_n, thtgrid, v_out, vm_out, vf_out, itheta_out):
    # this assumes block is for the same a and theta
    ia, ie, it = cuda.grid(3)
    
    #v_in_store  = cuda.shared.array((500,),f4)
    #vf_in_store = cuda.shared.array((500,),f4)
    #vm_in_store = cuda.shared.array((500,),f4)
    
    
    na = v_y.shape[0]
    ne = v_y.shape[1]
    nt = v_y.shape[2]
    
    
    if ia < na and ie < ne and it < nt:
        #v_in_store[it] = v_y[ia,ie,it]
        #vf_in_store[it] = vf_y[ia,ie,it]
        #vm_in_store[it] = vm_y[ia,ie,it]
        
        #cuda.syncthreads()
        
        vf_no = vf_n[ia,ie,it]
        vm_no = vm_n[ia,ie,it]
        
        
        
        #if vf_in_store[it] >= vf_no and vm_in_store[it] >= vm_no:
        if vf_y[ia,ie,it] >= vf_no and vm_y[ia,ie,it] >= vm_no:
            itheta_out[ia,ie,it] = it
            return
        
        #if vf_in_store[it] < vf_no and vm_in_store[it] < vm_no:
        if vf_y[ia,ie,it] < vf_no and vm_y[ia,ie,it] < vm_no:
            itheta_out[ia,ie,it] = -1
            tht = thtgrid[it]
            v_out[ia,ie,it] = tht*vf_no + (1-tht)*vm_no
            vf_out[ia,ie,it] = vf_no
            vm_out[ia,ie,it] = vm_no
            return
        
        
#        
#        it_ren = -1
#        
#        found_increase = False
#        found_decrease = False
#        
#        
#        for it_inc in range(it+1,nt):
#            if (vf_in_store[it_inc] >= vf_no and vm_in_store[it_inc] >= vm_no):
#                found_increase = True
#                break
#        
#        for it_dec in range(it-1,-1,-1):
#            if (vf_in_store[it_dec] >= vf_no and vm_in_store[it_dec] >= vm_no):
#                found_decrease = True
#                break
#            
#        
#        if found_increase and found_decrease:
#            dist_increase = it_inc - it
#            dist_decrease = it - it_dec
#            
#            if dist_increase != dist_decrease:
#                it_ren = it_inc if dist_increase < dist_decrease else it_dec
#            else:
#                # tie breaker
#                # numba-cuda does not do abs so we do these dumb things
#                dist_mid_inc = it_inc - (nt/2)                
#                if dist_mid_inc < 0: dist_mid_inc = -dist_mid_inc
#                dist_mid_dec = it_dec - (nt/2)
#                if dist_mid_dec < 0: dist_mid_dec = -dist_mid_dec
#            
#        elif found_increase and not found_decrease:
#            it_ren = it_inc
#        elif found_decrease and not found_increase:
#            it_ren = it_dec
#        else:
#            it_ren = -1 # check this!
#            
#        # finally fill the values    
#            
#        if it_ren == -1:
#            tht = thtgrid[it]
#            v_out[ia,ie,it] = tht*vf_no + (1-tht)*vm_no
#            vf_out[ia,ie,it] = vf_no
#            vm_out[ia,ie,it] = vm_no
#            itheta_out[ia,ie,it] = -1
#        else:
#            
#            # rescaling
#            tht_old = thtgrid[it]
#            tht_new = thtgrid[it_ren]
#            factor = (1-tht_old)/(1-tht_new) if tht_old < tht_new else tht_old/tht_new
#            
#            v_out[ia,ie,it] = factor*v_in_store[it_ren]
#            vf_out[ia,ie,it] = vf_in_store[it_ren]
#            vm_out[ia,ie,it] = vm_in_store[it_ren]
#            itheta_out[ia,ie,it] = it_ren
            
        
        
        
        
        
        
    
    
    
    
    
                
                