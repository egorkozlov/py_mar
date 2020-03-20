#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 13:47:08 2020

@author: egorkozlov
"""
import numpy as np
from aux_routines import first_true, last_true
from numba import njit, vectorize
from gridvec import VecOnGrid


def v_ren_new(setup,V,marriage,t,return_extra=False,return_vdiv_only=False,rescale=True):
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
    
    
    ind, izf, izm, ipsi = setup.all_indices(t+1)
    
    zfgrid = setup.exo_grids['Female, single'][t+1]
    zmgrid = setup.exo_grids['Male, single'][t+1]
    
    share=(np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) )
    relat=np.ones(share.shape)*0.5
    income_share_f=(1.0*share+0.0*relat).squeeze()
    #income_share_f = (np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) ).squeeze()
    
    share_f, share_m = dc.shares_if_split(income_share_f)
   
    
    
    
    # this is the part for bilateral divorce
    a_fem, a_mal = share_f[None,:]*setup.agrid_c[:,None], share_m[None,:]*setup.agrid_c[:,None]
    aleft_c = a_fem + a_mal
    iadiv_fem, iadiv_mal = [np.minimum(np.searchsorted(setup.agrid_s,x),setup.na-1)
                                        for x in [a_fem, a_mal]]
    na_s, nexo, ntheta = setup.na, ind.size, setup.ntheta_fine
    iadiv_fem_full, iadiv_mal_full = [np.broadcast_to(x[...,None],(na_s,nexo,ntheta)) 
                                        for x in [iadiv_fem, iadiv_mal]]
    aleft_c_full = np.broadcast_to(aleft_c[...,None],(na_s,nexo,ntheta))
    
    vf_all_s = V['Female, single']['V'][:,izf]
    vm_all_s = V['Male, single']['V'][:,izm]
        
    
    
    
    sc = setup.agrid_c
    
    # values of divorce
    vf_n, vm_n = v_div_byshare(
        setup, dc, t, sc, share_f, share_m,
        V['Male, single']['V'], V['Female, single']['V'],
        izf, izm, cost_fem=dc.money_lost_f, cost_mal=dc.money_lost_m)
    
    
    
    
    if return_vdiv_only:
        return {'Value of Divorce, male': vm_n,
                'Value of Divorce, female': vf_n}
    
    
    assert vf_n.ndim == vm_n.ndim == 2
    
    
    
    
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
        switch = (v_y_mar> v_y_coh)
        
        v_y = switch*v_y_mar + (~switch)*v_y_coh
        vf_y = switch*vf_y_mar + (~switch)*vf_y_coh
        vm_y = switch*vm_y_mar + (~switch)*vm_y_coh
        
    vf_n, vm_n = [x.astype(v_y.dtype) for x in (vf_n,vm_n)] # type conversion
    
    result  = v_ren_core_interp(setup,v_y, vf_y, vm_y, vf_n, vm_n, is_unil, rescale=rescale)
    
    result2 = ren_loop_wrap(setup, v_y, vf_y, vm_y, vf_n, vm_n, is_unil, rescale=rescale)
    
    
    from renegotiation_bilateral import ren_bilateral_wrap
    
    result3 = ren_bilateral_wrap(setup,v_y,vf_y,vm_y,vf_n,vm_n,vf_all_s,vm_all_s,aleft_c_full,                       
                       iadiv_fem_full,iadiv_mal_full,rescale=True)
    
    
    
    assert all([np.allclose(x,y) for x,y in zip(result['Values'],result2['Values'])])
    
    
    if not marriage:
        result['Cohabitation preferred to Marriage'] = ~switch
        
        
    
        
    extra = {'Values':result['Values'],
             'Value of Divorce, male': vm_n, 'Value of Divorce, female': vf_n}
    
    if not return_extra:
        return result
    else:
        return result, extra
    
    


def v_no_ren(setup,V,marriage,t):
    
    # this works live v_ren_new but does not actually run renegotiation
    
    expnd = lambda x : setup.v_thetagrid_fine.apply(x,axis=2)
    
    if marriage:
        # if couple is married already
        v_y = expnd(V['Couple, M']['V'])
        vf_y = expnd(V['Couple, M']['VF'])
        vm_y = expnd(V['Couple, M']['VM'])
    else:
        # stay in cohabitation
        v_y = expnd(V['Couple, C']['V'])
        vf_y = expnd(V['Couple, C']['VF'])
        vm_y = expnd(V['Couple, C']['VM'])
        switch = np.full(vf_y.shape,False,dtype=np.bool)
     
        
    def r(x): return x.astype(np.float32)
    
    shape_notheta = v_y.shape[:-1]
    yes = np.full(shape_notheta,True)
    ntheta = setup.ntheta_fine
    i_theta_out = np.broadcast_to(np.arange(ntheta,dtype=np.int16)[None,None,:],v_y.shape).copy()
        
    vf_n, vm_n = np.full((2,) + shape_notheta,-np.inf,dtype=np.float32)
    
    result =  {'Decision': yes, 'thetas': i_theta_out,
            'Values': (r(v_y), r(vf_y), r(vm_y)),'Divorce':(vf_n,vm_n)}
    
    if not marriage:
        result['Cohabitation preferred to Marriage'] = ~switch
        
        
        
    return result



def v_div_byshare(setup,dc,t,sc,share_fem,share_mal,Vmale,Vfemale,izf,izm,cost_fem=0.0,cost_mal=0.0):
    # this produces value of divorce for gridpoints given possibly different
    # shares of how assets are divided. 
    # Returns Vf_divorce, Vm_divorce -- values of singles in case of divorce
    # matched to the gridpionts for couples
    
    # optional cost_fem and cost_mal are monetary costs of divorce
    
    
    shrs = [0.2,0.35,0.5,0.65,0.8]  # grid on possible assets divisions    
    shp  =  (sc.size,izm.size,len(shrs))  
    Vm_divorce_M = np.zeros(shp) 
    Vf_divorce_M = np.zeros(shp)
    
    # find utilities of divorce for different divisions of assets
    for i, shr in enumerate(shrs):
        sv_m = VecOnGrid(setup.agrid_s, shr*sc - cost_mal)
        sv_f = sv_m if cost_fem == cost_mal else VecOnGrid(setup.agrid_s,shr*sc - cost_fem)
        
        Vm_divorce_M[...,i] = sv_m.apply(Vmale,    axis=0,take=(1,izm),reshape_i=True) - dc.u_lost_m
        Vf_divorce_M[...,i] = sv_f.apply(Vfemale,  axis=0,take=(1,izf),reshape_i=True) - dc.u_lost_f
    
    # share of assets that goes to the female
    # this has many repetative values but it turns out it does not matter much
    
    fem_gets = VecOnGrid(np.array(shrs),share_fem)
    mal_gets = VecOnGrid(np.array(shrs),share_mal)
    
    i_fem = fem_gets.i
    wn_fem = fem_gets.wnext
    
    i_mal = mal_gets.i
    wn_mal = mal_gets.wnext
    
    inds_exo = np.arange(setup.pars['nexo_t'][t+1])
    
    
    
    Vf_divorce = (1-wn_fem[None,:])*Vf_divorce_M[:,inds_exo,i_fem] + \
                wn_fem[None,:]*Vf_divorce_M[:,inds_exo,i_fem+1]
    
    Vm_divorce = (1-wn_mal[None,:])*Vm_divorce_M[:,inds_exo,i_mal] + \
                wn_mal[None,:]*Vm_divorce_M[:,inds_exo,i_mal+1]
                
    
                
    return Vf_divorce, Vm_divorce





def v_ren_core_interp(setup,vy,vfy,vmy,vf_n,vm_n,unilateral,show_sc=False,rescale=False):
    # this takes values of value functions (interpolated on fine grid)
    # and does discrete 
    # version of renegotiation.
    
    
    # compute the surplus
    
    
    
    sf_expand = vfy - vf_n[...,None]
    sm_expand = vmy - vm_n[...,None]
    
    exp_shape = sf_expand.shape
    
    
    # compute couple's value of divroce
    # make large arrays with values for each theta
    
    tgrid = setup.thetagrid_fine[None,None,:]
    ntheta = tgrid.size
    vf_div_full = np.broadcast_to(vf_n[...,None],exp_shape)
    vm_div_full = np.broadcast_to(vm_n[...,None],exp_shape)
    v_div_full = vf_div_full*tgrid + vm_div_full*(1-tgrid)
    
    
    
    # now do divorce
    
    
    
    v_out, vf_out, vm_out = vy.copy(), vfy.copy(), vmy.copy()
    i_theta_out = np.broadcast_to(np.arange(ntheta,dtype=np.int16)[None,None,:],exp_shape).copy()
        
    
    if not unilateral:
        # simple procedure
        no = (sf_expand<0) & (sm_expand<0)
        yes = ~no
        v_out[no] = v_div_full[no]
        vf_out[no] = vf_div_full[no]
        vm_out[no] = vm_div_full[no] # this has full size (has theta on the last axis)
        # therefore no ,:
        
        i_theta_out[no] = -1
        
        def r(x): return x.astype(np.float32)
        
        return {'Decision': yes, 'thetas': i_theta_out,
                'Values': (r(v_out), r(vf_out), r(vm_out)),'Divorce':(vf_n,vm_n)}
        
    # the rest handles unilateral divorce
    
    
    
    # compute where each agent agrees
    i_sf_expand = (sf_expand >= 0)
    i_sm_expand = (sm_expand >= 0)
    
    
    # agreement for all theta
    agree = (i_sf_expand) & (i_sm_expand)
    # any agreement     
    yes = np.any(agree,axis=-1)
    
    
    
    # then we divide the agreement points for single crossing and 
    # non-single crossing, as we need to handle renegotiation differently
    #
    # check for single crossing
    # signle crossing from false to true
    d_sf = np.sum(np.diff(i_sf_expand.astype(int),axis=-1),axis=-1) 
    n_sf = np.sum(np.abs(np.diff(i_sf_expand.astype(int),axis=-1)),axis=-1) 
    sc_f = ((d_sf == 1) & (n_sf==1)) | ((d_sf == 0) & (n_sf==0))
    # single crossing from true to false
    d_sm = np.sum(np.diff(i_sm_expand.astype(int),axis=-1),axis=-1)
    n_sm = np.sum(np.abs(np.diff(i_sm_expand.astype(int),axis=-1)),axis=-1)
    sc_m = ((d_sm == -1) & (n_sm == 1)) | ((d_sm == 0) & (n_sm==0))
    sc = (sc_f) & (sc_m)
    
   
    
    # agreement + single crossing
    yes_sc = (yes) & (sc)
    # agreement +  non-single crossing:
    yes_nsc = (yes) & ~(sc)         
    # disagreement
    no = ~(yes)
    
    share_sc = np.mean(yes_sc)
    share_nsc = np.mean(yes_nsc)
    
    if (share_nsc > 0) & show_sc: print('Not single crossing in {}, singe crossing in {} cases'.format(share_nsc,share_sc))
    
    # disagreement values
    v_out[no,:]  = v_div_full[no,:]
    vf_out[no,:] = vf_div_full[no,:]
    vm_out[no,:] = vm_div_full[no,:]
    i_theta_out[no,:] = -1
    
    # agreement values
    for yes_i, solver_i in zip([yes_sc,yes_nsc],[ind_sc,ind_no_sc]):
        if not np.any(yes_i): continue
            
        agree_this = agree[yes_i,:] # this is large matrix (??? x ntheta)
                                         # where (???) is all combinations of 
                                         # couple's characteristics that 
                                         # give agreement + single crossing
                                         
        inds = solver_i(agree_this)  # finds indices of nearest positive element
                                    # on a fine grid for theta
                                    # this is the most substantial part
        
        i_theta_out[yes_i,:] = inds # indexing is a bit mad :(
        
        v_out[yes_i,:] = np.take_along_axis(vy[yes_i,:],inds,axis=1)
        vf_out[yes_i,:] = np.take_along_axis(vfy[yes_i,:],inds,axis=1)
        vm_out[yes_i,:] = np.take_along_axis(vmy[yes_i,:],inds,axis=1)
        
        
    if not np.all(vf_out>=vf_div_full - 1e-4):
        print('Warning: f is broken is {} cases'.format(np.sum(vf_out<=vf_div_full - 1e-4)))
        raise Exception('regenotiation problems...')
        
    if not np.all(vm_out>=vm_div_full - 1e-4):
        print('Warning: m is broken is {} cases'.format(np.sum(vm_out<=vm_div_full - 1e-4)))
        raise Exception('regenotiation problems...')
    
    def r(x): return x.astype(np.float32)
    
   
    
    if rescale:
        theta_orig = np.broadcast_to(tgrid,i_theta_out.shape)
        theta_new  = setup.thetagrid_fine[i_theta_out]
        theta_new[no,:] = theta_orig[no,:] # this fixed divorced values
        factor = np.maximum((1-theta_orig)/(1-theta_new), theta_orig / theta_new)
        v_out_resc = factor*v_out    
        assert np.all(factor>=1)
        assert np.allclose(v_out_resc[no,:],v_out[no,:])
        v_out = v_out_resc
        # TODO: do we need to rescale vf_out, vm_out?
    
    return {'Decision': yes, 'thetas': i_theta_out,
            'Values': (r(v_out), r(vf_out), r(vm_out)),'Divorce':(vf_n,vm_n)}


    
    
def ren_loop_wrap(setup,vy,vfy,vmy,vfn,vmn,is_unil,rescale=False):
    # v_ren_core_interp(setup,vy,vfy,vmy,vf_n,vm_n,unilateral,show_sc=False,rescale=False)
    tgrid = setup.thetagrid_fine
    
    assert is_unil
    
    vout, vfout, vmout, thetaout, yes, ithetaout = \
        ren_loop(vy,vfy,vmy,vfn,vmn,tgrid,rescale=rescale)
    
    def r(x): return x.astype(np.float32)        
    
    
    return {'Decision': yes, 'thetas': ithetaout,
            'Values': (r(vout), r(vfout), r(vmout)),'Divorce':(vfn,vmn)}
    

@njit
def ren_loop(vy,vfy,vmy,vfn,vmn,thtgrid,rescale=False):
    print('hi!')


    vfn = vfn.reshape(vfn.shape+(1,))
    vmn = vmn.reshape(vmn.shape+(1,))
    
    sf = vfy - vfn
    sm = vmy - vmn
    
    na, nexo, nt = vy.shape
    
    vout = vy.copy()
    vfout = vfy.copy()
    vmout = vmy.copy()
    
    yes = np.zeros((na,nexo),dtype=np.bool_)
    
    
    thetaout = -1*np.ones(vout.shape,dtype=np.float32)
    ithetaout = -1*np.ones(vout.shape,dtype=np.int16)
    
    
    for ia in range(na):
        for iexo in range(nexo):
            sf_i = sf[ia,iexo,:]
            sm_i = sm[ia,iexo,:]
            
            both = (sf_i >= 0) & (sm_i >= 0)
            
            
            if not np.any(both):
                # divorce
                vfout[ia,iexo,:] = vfn[ia,iexo,0]
                vmout[ia,iexo,:] = vmn[ia,iexo,0]
                for itheta in range(nt):
                    th = thtgrid[itheta]
                    vout[ia,iexo,itheta] = th*vfn[ia,iexo,0] + (1-th)*vmn[ia,iexo,0]
            
            else:
                # renegotiate
                yes[ia,iexo] = True
                numbers = np.nonzero(both)[0]
                
                for itheta in range(nt):
                    if both[itheta]:
                        # status quo
                        thetaout[ia,iexo,itheta] = thtgrid[itheta]
                        ithetaout[ia,iexo,itheta] = itheta
                        continue
                    # if not both
                    in_closest = np.argmin(np.abs(numbers - itheta))
                    i_closest = numbers[in_closest]
                    ithetaout[ia,iexo,itheta] = i_closest
                    thetaout[ia,iexo,itheta] = thtgrid[i_closest]
                    
                    
                    if rescale:
                        tht_new = thtgrid[i_closest]
                        tht_old = thtgrid[itheta]
                        factor = np.maximum( (1-tht_old)/(1-tht_new), tht_old/tht_new )
                    else:
                        factor = 1
                    vout[ia,iexo,itheta] = factor*vy[ia,iexo,i_closest]
                    vfout[ia,iexo,itheta] = vfy[ia,iexo,i_closest] # no rescaling?
                    vmout[ia,iexo,itheta] = vmy[ia,iexo,i_closest]
                    
    if not np.all(vfout>=vfn - 1e-4):
        raise Exception('regenotiation problems...')
        
    if not np.all(vmout>=vmn - 1e-4):
        raise Exception('regenotiation problems...')                
    
    return vout, vfout, vmout, thetaout, yes, ithetaout




def ind_sc(i_pos):
    
        
        n_f = first_true(i_pos,axis=1)
        n_m =  last_true(i_pos,axis=1)
        
        npoints, nt = i_pos.shape
        
        assert not np.any(n_f==-1)
        assert not np.any(n_m==-1)
        
        inds = np.repeat(np.expand_dims(np.arange(nt),0),npoints,axis=0) # repeat
        n_f_bc = np.broadcast_to(n_f[:,None],(npoints, nt))
        n_m_bc = np.broadcast_to(n_m[:,None],(npoints, nt))
        
        i_f_less = (np.arange(nt)[None,:] < n_f[:,None])
        inds[i_f_less] = n_f_bc[i_f_less]
        
        i_m_more = (np.arange(nt)[None,:] > n_m[:,None])
        inds[i_m_more] = n_m_bc[i_m_more]
        
        return inds
    
    
@njit
def ind_no_sc(i_pos):
    # this uses pretty obscure loop
    nrows, ncol = i_pos.shape
    aran = np.arange(ncol)
    inds = np.empty((nrows,ncol),dtype=np.int32)
    
    for i in range(nrows):
        inds[i,:] = aran
        
        i_agree = i_pos[i,:]
        
        where_inds = np.nonzero(i_agree)[0]
        
        i_first_true = where_inds.min()
        i_last_true = where_inds.max()
        
        i_right = i_last_true
        
        j = 0
        while j < ncol:
            if i_agree[j]: 
                i_right = i_last_true
                j += 1
                continue
            if j < i_first_true:
                inds[i,j] = i_first_true
                j += 1
                continue
            if j > i_last_true:
                inds[i,j] = i_last_true
                j += 1
                continue
            
            # we are past region with positive values
            # and there is some region to the right with positive values
            
            
            # go right
            for jp in range(j+1,i_right+1,1):
                if i_agree[jp]:
                    i_right = jp
                    break
            
            # we found index of right
            
            i_left = j-1
            
            assert i_right > i_left
            
            for j in range(j,i_right):
                if j - i_left <= i_right - j:
                    inds[i,j] = i_left
                else:
                    inds[i,j] = i_right
            
            j+=1
            
            # we are here if above i_first_true and below i_last_true and 
            # i_agree[j] is false
                
            
            
    return inds
