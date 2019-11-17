#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for renegotiation and marriage

Possibly inefficient but very scalable is the goal

"""


#from trans_unif import transition_uniform
import numpy as np
from aux_routines import first_true, last_true
from numba import njit
from gridvec import VecOnGrid


# these are new routines for renegotiation

def v_ren_new(setup,V,marriage,t,sc=None,ind_or_inds=None,interpolate=True,combine=True,return_all=False):
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
        des='Couple, M'
    else:
        dc = setup.sep_costs
        is_unil = dc.unilateral_divorce # whether to do unilateral divorce at all
        des='Couple, C'
    
    
    if sc is None:
        sc = setup.agrid_c # savings of couple are given by the agrid
    
    if ind_or_inds is not None:
        ind, izf, izm, ipsi = setup.all_indices(ind_or_inds)
    else:
        ind, izf, izm, ipsi = setup.all_indices()
    
    
    assert (dc.money_lost_f_ez == 0 and dc.money_lost_m_ez == 0), 'not implemented yet'

    
    zfgrid = setup.exo_grids['Female, single'][t]
    zmgrid = setup.exo_grids['Male, single'][t]
    
    
    
    income_share_f = (np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) ).squeeze()
    
    share_f, share_m = dc.shares_if_split(income_share_f)
    
    
    
    if combine:
        # call separate routine for finding divorce outcomes depending on grid
        # position
        Vf_divorce, Vm_divorce = v_div_byshare(
            setup, dc, t, sc, share_f, share_m,
            V['Male, single']['V'], V['Female, single']['V'],
            izf, izm, cost_fem=dc.money_lost_f, cost_mal=dc.money_lost_m)
        
        assert Vf_divorce.ndim == Vm_divorce.ndim == 2
        
    else:
        # literally compute how assets are divided and apply
        # everything is 1-dimensional so no need for sophisticated things
        sm_v = VecOnGrid(setup.agrid_s,share_m*sc - dc.money_lost_m ,trim=True)
        sf_v = VecOnGrid(setup.agrid_s,share_f*sc - dc.money_lost_f ,trim=True)
        
        Vm_divorce = sm_v.apply(V['Male, single']['V'],  axis=0,take=(1,izm),reshape_i=combine)- dc.u_lost_m
        Vf_divorce = sf_v.apply(V['Female, single']['V'],axis=0,take=(1,izf),reshape_i=combine)- dc.u_lost_f
        
        try:
            assert Vm_divorce.ndim == Vf_divorce.ndim == 1
        except:
            print((Vm_divorce.ndim, Vm_divorce, Vf_divorce.ndim, Vf_divorce))
            assert False
    
    
    
    
    sc_v = VecOnGrid(setup.agrid_c,sc,trim=True)

    Vval_postren, VMval_postren, VFval_postren = \
        (sc_v.apply(v,axis=0,take=(1,ind),reshape_i=combine)
            for v in (V[des]['V'], V[des]['VM'], V[des]['VF']))
    
    try:
        assert Vval_postren.shape[:-1] == Vm_divorce.shape
    except:
        print(Vval_postren.shape,Vm_divorce.shape)
        assert False
    
    
    result  = v_ren_core_interp(setup,Vval_postren, VFval_postren, VMval_postren, Vf_divorce, Vm_divorce)
    
    
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
        sv_f = sv_m if cost_fem == cost_mal else VecOnGrid(setup.agrid,shr*sc - cost_fem)
        
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
    
    inds_exo = np.arange(setup.pars['nexo'])
    
    Vf_divorce = (1-wn_fem[None,:])*Vf_divorce_M[:,inds_exo,i_fem] + \
                wn_fem[None,:]*Vf_divorce_M[:,inds_exo,i_fem+1]
    
    Vm_divorce = (1-wn_mal[None,:])*Vm_divorce_M[:,inds_exo,i_mal] + \
                wn_mal[None,:]*Vm_divorce_M[:,inds_exo,i_mal+1]
                
    
                
    return Vf_divorce, Vm_divorce




def v_ren_core_interp(setup,v_y,vf_y,vm_y,vf_n,vm_n):
    # compute the surplus
    sf = vf_y - vf_n[...,None]
    sm = vm_y - vm_n[...,None]
    
    # interpolate it
    tgf = setup.v_thetagrid_fine
    sf_expand = tgf.apply(sf,axis=sf.ndim-1)
    sm_expand = tgf.apply(sm,axis=sm.ndim-1)
    
    
    
    
    
    
    
    # compute where each agent agrees
    i_sf_expand = (sf_expand >= 0)
    i_sm_expand = (sm_expand >= 0)
    
    
    # check for single crossing
    # signle crossing from false to true
    d_sf = np.sum(np.diff(i_sf_expand.astype(int),axis=-1),axis=-1) 
    sc_f = (d_sf == 1) | (d_sf == 0)
    # single crossing from true to false
    d_sm = np.sum(np.diff(i_sm_expand.astype(int),axis=-1),axis=-1)
    sc_m = (d_sm == -1) | (d_sm == 0)
    sc = (sc_f) & (sc_m)
    
    
    # agreement for all theta
    agree = (i_sf_expand) & (i_sm_expand)
    # any agreement     
    yes = np.any(agree,axis=-1)
    
    # then we split this by three regions: agreement + single crossing
    yes_sc = (yes) & (sc)
    # agreement +  non-single crossing:
    yes_nsc = (yes) & ~(sc)   
    # disagreement
    no = ~(yes)
    # the reason is that single crossing has much simpler algorithm
    
    # these things still have nice shape
    
    
    full_shape = sf.shape
    
    
    # compute couple's value of divroce
    # make large arrays with values for each theta
    tgrid = setup.thetagrid.reshape((1,)*(len(full_shape)-1) + (setup.ntheta,))
    vf_div_full = np.broadcast_to(vf_n[...,None],full_shape)
    vm_div_full = np.broadcast_to(vm_n[...,None],full_shape)
    v_div_full = vf_div_full*tgrid + vm_div_full*(1-tgrid)
    #theta_full = np.broadcast_to(tgrid,full_shape)
    
    
    v_out, vf_out, vm_out = v_y.copy(), vf_y.copy(), vm_y.copy()
    
    theta_out = -1*np.ones(full_shape,dtype=np.float32) # resulting thetas, can be off grid    
    i_theta_out = -1*np.ones(full_shape,dtype=np.int16) # grid position 
    wn_theta_out = np.zeros(full_shape,dtype=np.float32) # weight of the next point
    
    
    # fill disagreement
    v_out[no,:] = v_div_full[no,:]
    vf_out[no,:] = vf_div_full[no,:]
    vm_out[no,:] = vm_div_full[no,:]
    
    
    
    
    # renegotiation for single crossing points
    # note that this will be reshaped
    
    for yes_i, solver_i, solver_j in zip([yes_sc,yes_nsc],[ind_sc,ind_no_sc],[ind_no_sc,ind_sc]):
        if not np.any(yes_i): continue
            
        agree_this = agree[yes_i,:] # this is large matrix (??? x ntheta)
                                         # where (???) is all combinations of 
                                         # couple's characteristics that 
                                         # give agreement + single crossing
                                         
        #thetagrid_fine = setup.thetagrid_fine
        inds = solver_i(agree_this)  # finds indices of nearest positive element
                                    # on a fine grid for theta
        inds2 = solver_j(agree_this)
        assert np.allclose(inds,inds2)
        # then converts them to indices on regular grid
        theta_result = setup.v_thetagrid_fine.val[inds][:,setup.theta_orig_on_fine]
        i = setup.v_thetagrid_fine.i[inds][:,setup.theta_orig_on_fine]
        wn = setup.v_thetagrid_fine.wnext[inds][:,setup.theta_orig_on_fine]
        # and write them to the objects we care about
        
        theta_out[yes_i,:] = theta_result
        i_theta_out[yes_i,:] = i # indexing is a bit mad :(
        wn_theta_out[yes_i,:] = wn 
        
        
        fout = lambda x : (1-wn)*np.take_along_axis(x[yes_i,:],i,1) + \
                            wn*np.take_along_axis(x[yes_i,:],i+1,1)
        
        
        v_out[yes_i,:] = fout(v_y)
        vf_out[yes_i,:] = fout(vf_y)
        vm_out[yes_i,:] = fout(vm_y)
        
        
        
    assert np.all(vf_out>=vf_div_full - 1e-4)
    assert np.all(vm_out>=vm_div_full - 1e-4)
    assert np.all(theta_out[yes,:] > 0)
    assert not np.any(yes_nsc), 'Single crossing does not hold!' # FIXME: remove this later
    
    
    
    return {'Decision': yes, 'thetas': (theta_out, i_theta_out, wn_theta_out),
                'Values': (v_out, vf_out,vm_out)}




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




