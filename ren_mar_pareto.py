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
    inde=np.linspace(0,setup.ntheta-1,setup.ntheta,dtype=np.int16) 
     
    zfgrid = setup.exo_grids['Female, single'][t+1] 
    zmgrid = setup.exo_grids['Male, single'][t+1] 
     
    #share=(np.exp(zfgrid[izf]+setup.pars['f_wage_trend'][t]) / ( np.exp(zmgrid[izm]+setup.pars['m_wage_trend'][t]) + np.exp(zfgrid[izf]+setup.pars['f_wage_trend'][t]) ) ) 
    #relat=np.ones(share.shape)*0.5 
    #income_share_f=(1.0*share+0.0*relat).squeeze() 
    #income_share_f = (np.exp(zfgrid[izf]) / ( np.exp(zmgrid[izm]) + np.exp(zfgrid[izf]) ) ).squeeze() 
     
     
    #Get the share of female assets upon separation/divorce 
    share_f= dc.shares_if_split_theta(setup,setup.thetagrid)#dc.shares_if_split(income_share_f) 

    sc = setup.agrid_c 
     
    # values of divorce 
    vf_n, vm_n = v_div_byshare( 
        setup, dc, t, sc, share_f,inde, 
        V['Male, single']['V'], V['Female, single']['V'], 
        izf, izm, cost_fem=dc.money_lost_f, cost_mal=dc.money_lost_m) 
     
     
    if return_vdiv_only: 
        return {'Value of Divorce, male': vm_n, 
                'Value of Divorce, female': vf_n} 
     
     
    assert vf_n.ndim == vm_n.ndim == 3 
     
     
     
     
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
     
 
    return {'Values': (vfout, vmout), 'NBS': nbsout, 'theta': ithetaout, 'Decision':agree} 
 
 
 
def v_no_mar(setup,t,V,icouple,ind_or_inds,*,female,marriage): 
    # emulates v_mar_igrid but with no marriage 
     
     
    ind, izf, izm, ipsi = setup.all_indices(t,ind_or_inds) 
     
    vmout, vfout = V['Male, single']['V'][:,izm], V['Female, single']['V'][:,izf] 
     
     
    nbsout = np.zeros_like(vmout,dtype=np.float32) 
    ithetaout = -np.ones_like(vmout,dtype=np.int16) 
    agree = np.full_like(vmout,False,dtype=np.bool) 
     
    return {'Values': (vfout, vmout), 'NBS': nbsout, 'theta': ithetaout, 'Decision':agree} 
 
def v_div_byshare(setup,dc,t,sc,share_fem,inde,Vmale,Vfemale,izf,izm,cost_fem=0.0,cost_mal=0.0): 
    # this produces value of divorce for gridpoints given possibly different 
    # shares of how assets are divided.  
    # Returns Vf_divorce, Vm_divorce -- values of singles in case of divorce 
    # matched to the gridpionts for couples 
     
    # optional cost_fem and cost_mal are monetary costs of divorce 
     
     
    shrs = share_fem  # grid on possible assets divisions femle based     
    shp  =  (sc.size,izm.size,len(shrs))   
    Vm_divorce_M = np.zeros(shp)  
    Vf_divorce_M = np.zeros(shp) 
     
     
    Vmale_temp=np.zeros((len(setup.agrid_c),izm.size)) 
    Vfemale_temp=np.zeros((len(setup.agrid_c),izf.size)) 
    #Reshape first 
    for i in range(len(izf)):           
        Vfemale_temp[:,i] =Vfemale[:,izf[i]] 
         
    for i in range(len(izm)):           
        Vmale_temp[:,i] =Vmale[:,izm[i]] 
                 
                 
    # find utilities of divorce for different divisions of assets 
    for i, shr in enumerate(shrs): 
        sv_m = VecOnGrid(setup.agrid_s, (1-shr)*sc - cost_mal) 
        sv_f = VecOnGrid(setup.agrid_s,shr*sc - cost_fem) 
         
        wm=sv_m.wnext[...,None]*np.ones_like(izm) 
        wf=sv_f.wnext[...,None]*np.ones_like(izf) 
        Vm_divorce_M[:,:,i] = (1-wm)*Vmale_temp[sv_m.i,:] +(wm)*Vmale_temp[sv_m.i+1,:]- dc.u_lost_m 
        Vf_divorce_M[:,:,i] = (1-wf)*Vfemale_temp[sv_f.i,:] +(wf)*Vfemale_temp[sv_f.i+1,:] - dc.u_lost_f 
     
 
                 
    return Vf_divorce_M, Vm_divorce_M 
 
 
def v_ren_core_interp(setup,vy1,vfy1,vmy1,vf_n,vm_n,unilateral,show_sc=False,rescale=False): 
    # this takes values of value functions (interpolated on fine grid) 
    # and does discrete  
    # version of renegotiation. 
     
     
     # Expand to account that initial pareto weight matter 
    vfy=np.swapaxes(vfy1[...,None]*np.ones(setup.thetagrid.shape),2,3) 
    vmy=np.swapaxes(vmy1[...,None]*np.ones(setup.thetagrid.shape),2,3) 
    vy=np.swapaxes(vy1[...,None]*np.ones(setup.thetagrid.shape),2,3) 
     
     
    #setup.v_thetagrid_fine.apply(x,axis=2) 
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
     
         
        #The indexing below is a bit crazy but it seems to work 
        v_out1=v_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        vf_out1=vf_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        vm_out1=vm_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        i_theta_out1=i_theta_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        yes1=yes[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
 
         
        return {'Decision': yes1, 'thetas': i_theta_out1, 
                'Values': (r(v_out1), r(vf_out1), r(vm_out1)),'Divorce':(vf_n,vm_n)} 
         
    # the rest handles unilateral divorce 
     
     
     
    # compute where each agent agrees 
    i_sf_expand = (sf_expand >= 0) 
    i_sm_expand = (sm_expand >= 0) 
     
     
    # agreement for all theta 
    agree = (i_sf_expand) & (i_sm_expand) 
    # any agreement      
    yes = np.any(agree,axis=-1) 
    assert yes.ndim == 3 
     
     
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
         
         
        #The indexing below is a bit crazy but it seems to work 
        v_out1=v_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        vf_out1=vf_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        vm_out1=vm_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        i_theta_out1=i_theta_out[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1,range(setup.ntheta_fine)] 
        yes1=yes[:,:,setup.v_thetagrid_fine.i[range(setup.ntheta_fine)]+1] 
 
             
             
       
     
    return {'Decision': yes1, 'thetas': i_theta_out1, 
            'Values': (r(v_out1), r(vf_out1), r(vm_out1)),'Divorce':(vf_n,vm_n)} 
     
     
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
                     
'''
def ren_loop_bribe(vy,vfy,vmy,vfn,vmn,thtgrid,agrid_c,agrid_s): 
    vout = vy.copy() 
    vfout = vfy.copy() 
    vmout = vmy.copy() 
    
    na, nexo, nt = vy.shape 
     
    ithetaout = -1*np.ones(vout.shape,dtype=np.int32) 
    iassetsout_f = -2*np.ones(vout.shape,dtype=np.int32) 
    iassetsout_m = -2*np.ones(vout.shape,dtype=np.int32)
    
    ia_single_max_arr = np.searchsorted(agrid_c,agrid_s)

    
    for ia in range(na):
        for iexo in range(nexo):
            for it in range(nt):
                if vfy[ia,iexo,it] >= vfn[ia,iexo] and vmy[ia,iexo,it] >= vmn[ia,iexo]:
                    ithetaout[ia,iexo,it] = it
                    # do not fill iassets
                    continue # check if continue is valid
                
                if vfy[ia,iexo,it] < vfn[ia,iexo] and vmy[ia,iexo,it] < vmn[ia,iexo]:
                    # do not fill iassets
                    vfout[ia,iexo,it] = vfy[ia,iexo,it]
                    vmout[ia,iexo,it] = vmy[ia,iexo,it]
                    th = thtgrid[it]
                    vout[ia,iexo,it] = (1-th)*vmy[ia,iexo,it] + \
                                            th*vfy[ia,iexo,it]
                    continue # check if continue is valid
                
                if vfy[ia,iexo,it] < vfn[ia,iexo] or vmy[ia,iexo,it] < vmn[ia,iexo]:
                    # here only one constraint is binding
                    
                
'''                     
 
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
 
 
 
 
