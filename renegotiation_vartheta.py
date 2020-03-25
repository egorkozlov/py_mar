#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 19:33:42 2020

@author: egorkozlov
"""

import numpy as np
from aux_routines import first_true, last_true
from numba import njit, vectorize
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
        v_y  = expnd(V['Couple, M']['V'])
        vf_y = expnd(V['Couple, M']['VF'])
        vm_y = expnd(V['Couple, M']['VM'])
    else:
        # stay in cohabitation
        v_y_coh  = expnd(V['Couple, C']['V'])
        vf_y_coh = expnd(V['Couple, C']['VF'])
        vm_y_coh = expnd(V['Couple, C']['VM'])
        # switch to marriage
        v_y_mar  = expnd( V['Couple, M']['V'] )
        vf_y_mar = expnd( V['Couple, M']['VF'])
        vm_y_mar = expnd( V['Couple, M']['VM'])
        # switching criterion
        #switch = (vf_y_mar>vf_y_coh) & (vm_y_mar>vm_y_coh)
        switch = (v_y_mar >= v_y_coh)
        
        v_y = switch*v_y_mar + (~switch)*v_y_coh
        vf_y = switch*vf_y_mar + (~switch)*vf_y_coh
        vm_y = switch*vm_y_mar + (~switch)*vm_y_coh
        
    vf_n, vm_n = [x.astype(v_y.dtype) for x in (vf_n,vm_n)] # type conversion
    
    v_out, vf_out, vm_out, itheta_out = \
        v_ren_core(v_y, vf_y, vm_y, vf_n, vm_n, setup.thetagrid_fine, 
                   rescale=rescale)
        
    if marriage:
        itht = setup.v_thetagrid_fine.i
        wntht = setup.v_thetagrid_fine.wnext
        thtgrid = setup.thetagrid_fine
        v_out2, vf_out2, vm_out2, itheta_out2 = \
            v_ren_core_with_int(V['Couple, M']['V'],
                                V['Couple, M']['VF'], 
                                V['Couple, M']['VM'],
                                vf_n, vm_n,
                                itht, wntht, thtgrid, rescale = rescale)
        #try:
        assert np.all(itheta_out2==itheta_out)
        assert np.allclose(v_out2,v_out)
        assert np.allclose(vf_out2,vf_out)
        assert np.allclose(vm_out2,vm_out)
        print('worked!')
        #except:
        #    whr = np.where((itheta_out2!=itheta_out))
        #    print(np.where(whr))
        #    print(itheta_out[whr])
        #    print(itheta_out2[whr])
        #    print('failed')
        
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

                
                
@njit
def v_ren_core_with_int(v_y_ni, vf_y_ni, vm_y_ni, vf_n, vm_n, itht, wntht, thtgrid, rescale = True):
    # this takes values with no interpolation and interpolates inside
    
    
    na, ne, nt_coarse = v_y_ni.shape
    nt = thtgrid.size
    
    shp = (na,ne,nt)
    
    v_out = np.empty(shp,dtype=np.float32)
    vm_out = np.empty(shp,dtype=np.float32)
    vf_out = np.empty(shp,dtype=np.float32)
    
    itheta_out = np.full(v_out.shape,-1,dtype=np.int16)
    
    
    f1 = np.float32(1)
    
    
    for ia in range(na):
        for ie in range(ne):
            
            
            
            v_opt = np.empty((nt,),dtype=np.float32)
            vf_opt = np.empty((nt,),dtype=np.float32)
            vm_opt = np.empty((nt,),dtype=np.float32)
            
            # this part does all interpolations and maximization
            for it in range(nt):
                it_c = itht[it]
                it_cp = it_c+1
                wn_c = wntht[it]
                wt_c = f1 - wn_c
                
                def wsum(x):
                    return x[ia,ie,it_c]*wt_c + x[ia,ie,it_cp]*wn_c
                
                v_opt[it] = wsum(v_y_ni)
                vf_opt[it] = wsum(vf_y_ni)
                vm_opt[it] = wsum(vm_y_ni)
                
            
            # this part actually does renegotiation
            for it in range(nt):
                
                
                vf_y = vf_opt[it]                
                vm_y = vm_opt[it]
                v_y = v_opt[it]
                
                vf_no = vf_n[ia,ie,it]
                vm_no = vm_n[ia,ie,it]
                
                if vf_y >= vf_no and vm_y >= vm_no:
                    # no search just fill the value
                    itheta_out[ia,ie,it] = it    
                    vf_out[ia,ie,it] = vf_y
                    vm_out[ia,ie,it] = vm_y
                    v_out[ia,ie,it] = v_y
                    continue
                    
                if vf_y < vf_no and vm_y < vm_no:
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
                    if (vf_opt[it_increase] >= vf_no and vm_opt[it_increase] >= vm_no):
                        found_increase = True
                        break
                
                
                
                for it_decrease in range(it-1,-1,-1):
                    if (vf_opt[it_decrease] >= vf_no and vm_opt[it_decrease] >= vm_no):
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
                    
                    
                    vf_y = vf_opt[it_ren]              
                    vm_y = vm_opt[it_ren]
                    v_y  =  v_opt[it_ren]
                    
                    v_out[ia,ie,it] = factor*v_y
                    vf_out[ia,ie,it] = vf_y
                    vm_out[ia,ie,it] = vm_y
                    itheta_out[ia,ie,it] = it_ren
                
    
    return v_out, vf_out, vm_out, itheta_out


                
@njit
def v_ren_core_two_opts_with_int(v_y_ni_0, v_y_ni_1, vf_y_ni_0, vf_y_ni_1, vm_y_ni_0, vm_y_ni_1, vf_n, vm_n, itht, wntht, thtgrid, rescale = True):
    # this takes values with no interpolation and interpolates inside
    # this also makes a choice of mar / coh
    # choice is based on comparing v_y_ni_0 vs v_y_ni_1 in the interpolated pt
    
    na, ne, nt_coarse = v_y_ni_0.shape
    nt = thtgrid.size
    
    shp = (na,ne,nt)
    
    v_out = np.empty(shp,dtype=np.float32)
    vm_out = np.empty(shp,dtype=np.float32)
    vf_out = np.empty(shp,dtype=np.float32)
    
    itheta_out = np.full(v_out.shape,-1,dtype=np.int16)
    ichoice_out = np.zeros(v_out.shape,dtype=np.bool)
    
    
    f1 = np.float32(1)
    
    
    for ia in range(na):
        for ie in range(ne):
            # first we form value functions and choices
            # then we do renegotiation
            # this saves lots of operations
            
            v_opt = np.empty((nt,),dtype=np.float32)
            vf_opt = np.empty((nt,),dtype=np.float32)
            vm_opt = np.empty((nt,),dtype=np.float32)
            
            # this part does all interpolations and maximization
            for it in range(nt):
                it_c = itht[it]
                it_cp = it_c+1
                wn_c = wntht[it]
                wt_c = f1 - wn_c
                
                def wsum(x):
                    return x[ia,ie,it_c]*wt_c + x[ia,ie,it_cp]*wn_c
                
                v_y_0 = wsum(v_y_ni_1)
                v_y_1 = wsum(v_y_ni_0)                
                pick_1 = (v_y_1 > v_y_0)
                
                if pick_1:
                    vf_opt[it] = wsum(vf_y_ni_0)
                    vm_opt[it] = wsum(vm_y_ni_0)
                    v_opt[it] = v_y_0
                else:
                    vf_opt[it] = wsum(vf_y_ni_1)
                    vm_opt[it] = wsum(vm_y_ni_1)
                    v_opt[it] = v_y_1
                
                ichoice_out[ia,ie,it] = pick_1
                
                
            
            for it in range(nt):
                
                
                vf_y = vf_opt[it]
                vm_y = vm_opt[it]                
                v_y = v_opt[it]
                
                
                vf_no = vf_n[ia,ie,it]
                vm_no = vm_n[ia,ie,it]
                
                if vf_y >= vf_no and vm_y >= vm_no:
                    # no search just fill the value
                    itheta_out[ia,ie,it] = it    
                    vf_out[ia,ie,it] = vf_y
                    vm_out[ia,ie,it] = vm_y
                    v_out[ia,ie,it] = v_y
                    continue
                    
                if vf_y < vf_no and vm_y < vm_no:
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
                    if (vf_opt[it_increase] >= vf_no and vm_opt[it_increase] >= vm_no):
                        found_increase = True
                        break
                
                
                
                for it_decrease in range(it-1,-1,-1):
                    if (vf_opt[it_decrease] >= vf_no and vm_opt[it_decrease] >= vm_no):
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
                    
                        
                    v_out[ia,ie,it] = factor*v_opt[it_ren]
                    vf_out[ia,ie,it] = vf_opt[it_ren]
                    vm_out[ia,ie,it] = vm_opt[it_ren]
                    itheta_out[ia,ie,it] = it_ren
                
    
    return v_out, vf_out, vm_out, itheta_out, ichoice_out

