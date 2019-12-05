#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for renegotiation and marriage

Possibly inefficient but very scalable is the goal

"""


#from trans_unif import transition_uniform
import numpy as np
from aux_routines import first_true, last_true
from numba import njit, vectorize
from gridvec import VecOnGrid


# these are new routines for renegotiation

##### main part

def v_ren_new(setup,V,marriage,t):
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
        des=['Couple, M']
    else:
        dc = setup.sep_costs
        is_unil = dc.unilateral_divorce # whether to do unilateral divorce at all
        des=['Couple, C','Couple, M']
    
    assert is_unil, 'only unilateral divorce is implemented'
    
    ind, izf, izm, ipsi = setup.all_indices(t)
    
    zfgrid = setup.exo_grids['Female, single'][t]
    zmgrid = setup.exo_grids['Male, single'][t]
    
    income_share_f = (np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) ).squeeze()
    
    share_f, share_m = dc.shares_if_split(income_share_f)
    
    sc = setup.agrid_c
    
    
    Vf_divorce, Vm_divorce = v_div_byshare(
        setup, dc, t, sc, share_f, share_m,
        V['Male, single']['V'], V['Female, single']['V'],
        izf, izm, cost_fem=dc.money_lost_f, cost_mal=dc.money_lost_m)
    
    assert Vf_divorce.ndim == Vm_divorce.ndim == 2
    
    Vval_postren,VMval_postren,VFval_postren=\
    np.zeros([3,len(des),np.size(V['Couple, M']['V'],0),np.size(V['Couple, M']['V'],1),np.size(V['Couple, M']['V'],2)])
    

        
    for ist,sname in enumerate(des):
  
        Vval_postren[ist,], VMval_postren[ist,], VFval_postren[ist,] = V[sname]['V'], V[sname]['VM'], V[sname]['VF']
    
    assert np.size(Vval_postren,1) == np.size(Vm_divorce,0)
    assert np.size(Vval_postren,2) == np.size(Vm_divorce,1)
    
    
    result  = v_ren_core_interp(setup,Vval_postren, VFval_postren, VMval_postren, Vf_divorce, Vm_divorce)
    
    
    return result




def v_mar_igrid(setup,t,V,icouple,ind_or_inds,*,female,marriage,interpolate=True,return_all=False):
    # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    # if return_all==False returns Vout_f, Vout_m, that are value functions
    # of male and female from entering this union
    # if return_all==True returns (Vout_f, Vout_m, ismar, thetaout, technical)
    # where ismar is marriage decision, thetaout is resulting theta and 
    # tuple technical contains less usable stuff (check v_newmar_core for it)
    #
    # combine = True creates matrix (n_s-by-n_inds)
    # combine = False assumed that n_s is the same shape as n_inds and creates
    # a flat array.
    
    
    if marriage:
        coup = 'Couple, M'
    else:
        coup = 'Couple, C'
    
    
    
    # import objects
    agrid_c = setup.agrid_c
    agrid_s = setup.agrid_s
    gamma = setup.pars['m_bargaining_weight']   
    
    
    VMval_single, VFval_single = V['Male, single']['V'], V['Female, single']['V']
    VMval_postren, VFval_postren = V[coup]['VM'][icouple,...], V[coup]['VF'][icouple,...]
    
    # type conversion
    #VMval_single, VFval_single, VMval_postren, VFval_postren = \
    #    [np.float32(x) for x in  
    #     (VMval_single, VFval_single, VMval_postren, VFval_postren)]
    
    
    
    # substantial part
    ind, izf, izm, ipsi = setup.all_indices(t,ind_or_inds)
    
    
    # using trim = True implicitly trims things on top
    # so if sf is 0.75*amax and sm is 0.75*amax then sc is 1*amax and not 1.5
    
    #sc = sf+sm # savings of couple
    s_partner = agrid_c[icouple] - agrid_s # we assume all points on grid
    
    
    # this implicitly trims negative or too large values
    s_partner_v = VecOnGrid(agrid_s,s_partner,trim=True) 
    
    
    # this applies them
    
    if female:
        Vfs = VFval_single[:,izf]
        Vms = s_partner_v.apply(VMval_single,axis=0,take=(1,izm))
    else:
        Vms = VMval_single[:,izm]
        Vfs = s_partner_v.apply(VFval_single,axis=0,take=(1,izf))
        
        
        
    
    
    expnd = lambda x : setup.v_thetagrid_fine.apply(x,axis=2)
    
    
    Vmm, Vfm = (expnd(x[:,ind,:]) for x in 
                     (VMval_postren,VFval_postren))
    
   
    ins = [Vfm,Vmm,Vfs,Vms,gamma]
    ins = [np.float32(x) for x in ins] # optional type conversion
    vfout, vmout, nbsout, agree, ithetaout = mar_mat(*ins)
    

    #assert np.allclose(vfout,vfout2)
    #assert np.allclose(vmout,vmout2)
    #assert np.allclose(nbsout,nbsout2)
    #assert np.allclose(ithetaout,ithetaout2)
    
    
    
    return {'Values': (vfout, vmout), 'NBS': nbsout, 'theta': ithetaout, 'Decision':agree}

##### technical part
    




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
    
    inds_exo = np.arange(setup.pars['nexo_t'][t])
    
    
    
    Vf_divorce = (1-wn_fem[None,:])*Vf_divorce_M[:,inds_exo,i_fem] + \
                wn_fem[None,:]*Vf_divorce_M[:,inds_exo,i_fem+1]
    
    Vm_divorce = (1-wn_mal[None,:])*Vm_divorce_M[:,inds_exo,i_mal] + \
                wn_mal[None,:]*Vm_divorce_M[:,inds_exo,i_mal+1]
                
    
                
    return Vf_divorce, Vm_divorce


def v_ren_core_interp(setup,v_y,vf_y,vm_y,vf_n,vm_n):
    # compute the surplus
    
    # interpolate it
    tgf = setup.v_thetagrid_fine
    
    #Get relevant variables if marriage or cohabitation
    v_y_check=np.copy(v_y)
    if np.size(v_y,0)==1:
        v_y=v_y[0]
        vf_y=vf_y[0]
        vm_y=vm_y[0] 
        
        
        vy, vfy, vmy = (tgf.apply(x,axis=2) for x in (v_y,vf_y,vm_y) )
    else:
        
        assert np.size(v_y,0)==2
        
        #TODO try when both agree
        
        #Cohabitation
        v_y0=v_y[0]
        vf_y0=vf_y[0]
        vm_y0=vm_y[0] 
        
        #Marriage
        v_y1=v_y[1]
        vf_y1=vf_y[1]
        vm_y1=vm_y[1] 
        
        
        vy0, vfy0, vmy0 = (tgf.apply(x,axis=2) for x in (v_y0,vf_y0,vm_y0) )
        vy1, vfy1, vmy1 = (tgf.apply(x,axis=2) for x in (v_y1,vf_y1,vm_y1) )
        
        #Foreach theta_fine, get marriage or cohabitation
        #inde = (vy0>=vy1)            
        indef = (vfy0<vfy1)
        indem = (vmy0<vmy1)
        inde=1-indem*indef
        vy=inde*vy0+(1-inde)*vy1
        vfy=inde*vfy0+(1-inde)*vfy1
        vmy=inde*vmy0+(1-inde)*vmy1
    
    
    
    sf_expand = vfy - vf_n[...,None]
    sm_expand = vmy - vm_n[...,None]
    
    exp_shape = sf_expand.shape
    
    
    
    
    
    # compute where each agent agrees
    i_sf_expand = (sf_expand >= 0)
    i_sm_expand = (sm_expand >= 0)
    
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
    
    
    share_sc = np.mean(yes_sc)
    share_nsc = np.mean(yes_nsc)
    
    if share_nsc > 0: print('Not single crossing in {}, singe crossing in {} cases'.format(share_nsc,share_sc))
    
    # these things still have nice shape
    
    # compute couple's value of divroce
    # make large arrays with values for each theta
    
    tgrid = tgf.val[None,None,:]
    vf_div_full = np.broadcast_to(vf_n[...,None],exp_shape)
    vm_div_full = np.broadcast_to(vm_n[...,None],exp_shape)
    v_div_full = vf_div_full*tgrid + vm_div_full*(1-tgrid)
    #theta_full = np.broadcast_to(tgrid,full_shape)
    
    
    v_out, vf_out, vm_out = vy.copy(), vfy.copy(), vmy.copy()
    
    i_theta_out = -1*np.ones(exp_shape,dtype=np.int16) # grid position 
    
    
    # fill disagreement
    
    v_out[no,:]  = v_div_full[no,:]
    vf_out[no,:] = vf_div_full[no,:]
    vm_out[no,:] = vm_div_full[no,:]
    
    
    
    
    # renegotiation for single crossing points
    # note that this will be reshaped
    
    for yes_i, solver_i in zip([yes_sc,yes_nsc],[ind_sc,ind_no_sc]):
        if not np.any(yes_i): continue
            
        agree_this = agree[yes_i,:] # this is large matrix (??? x ntheta)
                                         # where (???) is all combinations of 
                                         # couple's characteristics that 
                                         # give agreement + single crossing
                                         
        #thetagrid_fine = setup.thetagrid_fine
        inds = solver_i(agree_this)  # finds indices of nearest positive element
                                    # on a fine grid for theta
        #inds2 = solver_j(agree_this)
        #assert np.allclose(inds,inds2)
        
        i_theta_out[yes_i,:] = inds # indexing is a bit mad :(
        
        
        v_out[yes_i,:] = np.take_along_axis(vy[yes_i,:],inds,axis=1)
        vf_out[yes_i,:] = np.take_along_axis(vfy[yes_i,:],inds,axis=1)
        vm_out[yes_i,:] = np.take_along_axis(vmy[yes_i,:],inds,axis=1)
        
        
    try:
        assert np.all(vf_out>=vf_div_full - 1e-4)
    except:
        print('Warning: broken is {} cases'.format(np.sum(vf_out<=vf_div_full - 1e-4)))
        
    try:
        assert np.all(vm_out>=vm_div_full - 1e-4)
    except:
        print('Warning: broken is {} cases'.format(np.sum(vm_out<=vm_div_full - 1e-4)))
        
    
    #assert not np.any(yes_nsc), 'Single crossing does not hold!' # FIXME: remove this later
    
    if np.size(v_y_check,0)==1:
        return {'Decision': yes, 'thetas': i_theta_out,
                'Values': (v_out, vf_out, vm_out),'Divorce':(vf_div_full,vm_div_full)}
    else:
        return {'Decision': yes, 'thetas': i_theta_out,
                'Values': (v_out, vf_out, vm_out),'Cohabitation preferred to Marriage': inde}



@njit
def ren_loop(vy,vfy,vmy,vfn,vmn,thtgrid):
    print('hi!')

    sf = vfy - vfn
    sm = vmy - vmn
    
    na, nexo, nt = vy.shape
    
    vout = vy.copy()
    vfout = vfy.copy()
    vmout = vmy.copy()
    
    thetaout = -1*np.zeros(vout.shape,dtype=np.float32)
    
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
                
                numbers = np.nonzero(both)[0]
                
                for itheta in range(nt):
                    if both[itheta]:
                        # status quo
                        thetaout[ia,iexo,itheta] = thtgrid[itheta]
                        continue
                    # if not both
                    in_closest = np.argmin(np.abs(numbers - itheta))
                    i_closest = numbers[in_closest]
                    thetaout[ia,iexo,itheta] = thtgrid[i_closest]
                    vout[ia,iexo,itheta] = vy[ia,iexo,i_closest]
                    vfout[ia,iexo,itheta] = vfy[ia,iexo,i_closest]
                    vmout[ia,iexo,itheta] = vmy[ia,iexo,i_closest]
                    
    
    return vout, vfout, vmout, thetaout
                    
                    
                    
                    

@njit
def mar_loop(vfy,vmy,vfn,vmn,gamma):

    sf = vfy - np.expand_dims(vfn,vfn.ndim)
    sm = vmy - np.expand_dims(vmn,vmn.ndim)
    
    na, nexo, nt = vfy.shape
    
    #vout = np.zeros((vy.shape[:-1]),np.float32)
    vfout = vfn.copy()
    vmout = vmn.copy()
    nbsout = np.zeros(vfn.shape)
    
    ithetaout = -1*np.ones(vfout.shape,dtype=np.int32)
    agree = np.zeros(vfout.shape,dtype=np.bool)
    
    ntheta = vfy.shape[-1]
    
    for ia in range(na):
        for iexo in range(nexo):
            sf_i = sf[ia,iexo,:]
            sm_i = sm[ia,iexo,:]
            
            both = (sf_i > 0) & (sm_i > 0)
            
            good = np.any(both)
             
            agree[ia,iexo] = good
            
            
            if good:
                nbs = np.zeros(ntheta,dtype=np.float32)
                nbs[both] = (sf_i[both]**gamma) * (sm_i[both]**(1-gamma))
                i_best = nbs.argmax()
                nbs_best = nbs[i_best]
                assert nbs_best > 0
                ithetaout[ia,iexo] = i_best
                vfout[ia,iexo] = vfy[ia,iexo,i_best]
                vmout[ia,iexo] = vmy[ia,iexo,i_best]
                nbsout[ia,iexo] = nbs_best
                
    return vfout, vmout, nbsout, agree, ithetaout
            
            

@vectorize('float32(float32,float32,float32)')  
def nbs(x,y,gamma):
    if x > 0 and y > 0:
        return (x**gamma) * (y**(1-gamma))
    else:
        return 0
                        


def mar_mat(vfy,vmy,vfn,vmn,gamma):

    sf = vfy - np.expand_dims(vfn,vfn.ndim)
    sm = vmy - np.expand_dims(vmn,vmn.ndim)
    
    
    
    #vout = np.zeros((vy.shape[:-1]),np.float32)
    vfout = vfn.copy()
    vmout = vmn.copy()
    
    
    agree = (sf>0) & (sm>0)
    any_agree = np.any(agree,axis=-1)
    
    # this reshapes things
    n_agree = np.sum(any_agree)
    
    nbsout = np.zeros(vfn.shape,dtype=np.float32)
    ithetaout = -1*np.ones(vfn.shape,dtype=np.int32)
    
    
    if n_agree > 0:       
        
        sf_a = sf[any_agree,:]
        sm_a = sm[any_agree,:]
        nbs_a = np.zeros(sf_a.shape,dtype=np.float32)
        
        a_pos = (sf_a>0) & (sm_a>0)
        
        nbs_a[a_pos] = (sf_a[a_pos]**gamma) * (sm_a[a_pos]**(1-gamma))
        inds_best = np.argmax(nbs_a,axis=1)
        
        take = lambda x : np.take_along_axis(x,inds_best[:,None],axis=1).reshape((n_agree,))
        
        nbsout[any_agree] = take(nbs_a) 
        assert np.all(nbsout[any_agree] > 0)
        ithetaout[any_agree] = inds_best
        vfout[any_agree] = take(vfy[any_agree,:])
        vmout[any_agree] = take(vmy[any_agree,:])
        
        
        
        
    return vfout, vmout, nbsout, any_agree, ithetaout
        
        
        
    
    
    
    #nbsout = nbs(sf,sm,gamma)
    
    
    
    
    
    #ithetaout = -1*np.zeros(vfout.shape,dtype=np.float32)
    
    
    
                
    
    #return vfout, vmout, nbsout, ithetaout
    


 


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




