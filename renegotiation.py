#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects things relevant for renegotiation of couples
"""
from trans_unif import transition_uniform
import numpy as np
from aux_routines import first_true, last_true
from numba import jit


def v_last_period_renegotiated(setup,V,kappa=0.45,interpolate=False):
    # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    ism, psm = transition_uniform(setup.agrid,kappa*setup.agrid)
    isf, psf = transition_uniform(setup.agrid,kappa*setup.agrid)
    
    ind, izf, izm, ipsi = setup.all_indices()
    
    
    nt = setup.thetagrid.size
    
    VMval_single, VFval_single = V['Male, single']['V'], V['Female, single']['V']
    Vval_postren, VMval_postren, VFval_postren = V['Couple']['V'], V['Couple']['VM'], V['Couple']['VF']
    
    
    Vm_divorce = (VMval_single[ism,:]*psm[:,None] + VMval_single[ism+1,:]*(1-psm[:,None]))[:,izm,None]
    Vf_divorce = (VFval_single[isf,:]*psf[:,None] + VFval_single[isf+1,:]*(1-psf[:,None]))[:,izf,None]
    
    S_f = VFval_postren - Vf_divorce # surplus of female
    S_m = VMval_postren - Vm_divorce # surplus of male
    
    I_f = np.array(S_f > 0) # whether female agrees at this gridpoint
    I_m = np.array(S_m > 0) # whether male agrees at this gridpoint
    
    
    sq = (I_f & I_m)
    
    
    nf = first_true(I_f,axis=2)
    nm = last_true(I_m,axis=2)
    
    together = np.any(sq,axis=2)
    
    
    
    together_2 = (np.array(nf<=nm) & np.array(nf!=-1) & np.array(nm!=-1) )
    
    on_the_edge = (np.array(nf == nm+1) & np.array(nf!=-1) & np.array(nm!=-1) )[:,:,None]
    
    
    
    assert np.all(together_2 == together)
    
    
    NF = nf[:,:,None]
    NF_fix = np.array( NF == 0 ) # when first true is zero
    NF[NF_fix] = 1
    NF_ok = np.array( NF != -1 ) # when there is at least some true value
    NF_sc = (NF_ok & ~NF_fix) # whether it corresponds to sign change
    NF[~NF_ok] = nt-1
    
    
    
    NM = nm[:,:,None]
    NM_fix = np.array( NM == nt-1)
    NM[NM_fix] = nt-2
    NM_ok = np.array(NM != -1)
    NM[~NM_ok] = 0
    NM_sc = (NM_ok & ~NM_fix)
    
    
    # rf is theta grid poisition to the right of point where female surplus intersects zero
    # (so it is positive as surplus is increasing in theta)
    # lf is to the left (so it is negative)
    # rm is theta grid poisition to the right of point where male surplus intersects zero
    # (so it is negative)
    # lm is to the left (so it is positive)
    
    # if renegotiated female-induced: move to rf
    # if renegotiated male-induced: move to lm
    # lf and rm are technical for interpolation
    
    rf, lf, rm, lm  = NF, NF-1, NM+1, NM
    
    VF_rf, VF_lf, VF_rm, VF_lm  = (np.take_along_axis(VFval_postren,x,2) for x in (rf, lf, rm, lm))
    VM_rf, VM_lf, VM_rm, VM_lm =  (np.take_along_axis(VMval_postren,x,2) for x in (rf, lf, rm, lm))
    VC_rf, VC_lf, VC_rm, VC_lm =  (np.take_along_axis( Vval_postren,x,2) for x in (rf, lf, rm, lm))
    

    A_l, A_r = VF_lf - Vf_divorce, VF_rf - Vf_divorce
    B_l, B_r = VM_lm - Vm_divorce, VM_rm - Vm_divorce
    
    assert np.all(A_r[NF_sc]>0)
    assert np.all(A_l[NF_sc]<0)
    assert np.all(B_l[NM_sc]>0)
    assert np.all(B_r[NM_sc]<0)
    
    
    kr_f = -A_l/(A_r - A_l)    
    kr_m = B_l/(B_l - B_r) 
    kr_f[NF_fix] = 0.0
    kr_m[NF_fix] = 1.0
    
    
    # out of NF_sc these values do not make much sense
    
    assert np.max(np.abs( (A_l*(1-kr_f) + A_r*kr_f)[NF_sc] )) < 1e-4
    assert np.max(np.abs( (B_l*(1-kr_m) + B_r*kr_m)[NM_sc] )) < 1e-4
    
    assert np.all(kr_f[NF_sc]>=0) and np.all(kr_f[NF_sc]<=1)
    assert np.all(kr_m[NM_sc]>=0) and np.all(kr_m[NM_sc]<=1)
    
    
    
    
    on_the_edge_y = np.full_like(on_the_edge,False)
    on_the_edge_y[(np.array(kr_f <= kr_m) & on_the_edge)] = True
    
    if np.any(on_the_edge_y):
        assert np.all( (B_l*(1-kr_f) + B_r*kr_f)[on_the_edge_y] >= 0 ) 
        assert np.all( (A_l*(1-kr_m) + A_r*kr_m)[on_the_edge_y] >= 0 )
        
    
    
    
    
    if interpolate:
        together = (together | on_the_edge_y.squeeze())
    
    
    
    f_ren = (~I_f & together[:,:,None]) # female-initiated renegotiation
    m_ren = (~I_m & together[:,:,None]) # male-initiated renegotiation
     
    
    
    # I create new for preren
    VF_out = np.copy(VFval_postren)
    VM_out = np.copy(VMval_postren)
    V_out  = np.copy(Vval_postren)
    
    if not interpolate:
        Vf_ifren_f = VF_rf
        Vf_ifren_m = VF_lm
        Vm_ifren_f = VM_rf
        Vm_ifren_m = VM_lm
        V_ifren_f  = VC_rf
        V_ifren_m  = VC_lm
    
    else:
        
        # when no agreement kr_f and kr_m are bullshit
        # this is alright as we do use only those that
        # are in f_ren and m_ren regions where there is an agreement
        
        Vf_ifren_f = (1-kr_f)*VF_lf + kr_f*VF_rf
        Vm_ifren_f = (1-kr_f)*VM_lf + kr_f*VM_rf
        V_ifren_f  = (1-kr_f)*VC_lf + kr_f*VC_rf
        
        Vf_ifren_m = (1-kr_m)*VF_lm + kr_m*VF_rm
        Vm_ifren_m = (1-kr_m)*VM_lm + kr_m*VM_rm
        V_ifren_m  = (1-kr_m)*VC_lm + kr_m*VC_rm
        
        
        
        
        
        assert np.all( np.abs( Vf_ifren_f - Vf_divorce)[NF_sc] < 1e-4 )
        assert np.all( np.abs( Vm_ifren_m - Vm_divorce)[NM_sc] < 1e-4 )
        
        
        
        
    bt = lambda x : np.broadcast_to(x, VF_out.shape) # mad skillz
    # this assumed VF_out and VM_out have the same shape
    
    t_stretch = bt(setup.thetagrid[None,None,:])
    
    
    
    bool_divorce = bt(~together[:,:,None])
    VF_out[bool_divorce] = bt(Vf_divorce)[bool_divorce]
    VM_out[bool_divorce] = bt(Vm_divorce)[bool_divorce]
    
    V_out[bool_divorce] = t_stretch[bool_divorce]*VF_out[bool_divorce] + \
                      (1-t_stretch[bool_divorce])*VM_out[bool_divorce] # check me please
    
    VF_out[f_ren] = bt(Vf_ifren_f)[f_ren]
    VM_out[f_ren] = bt(Vm_ifren_f)[f_ren]
    V_out[f_ren]  = bt( V_ifren_f)[f_ren]
    
    VF_out[m_ren] = bt(Vf_ifren_m)[m_ren]
    VM_out[m_ren] = bt(Vm_ifren_m)[m_ren]
    V_out[m_ren]  = bt( V_ifren_m)[m_ren]
    
    
    
    
    assert np.all(VF_out >= Vf_divorce - 1e-4)
    assert np.all(VM_out >= Vm_divorce - 1e-4)
    
    
    return V_out, VF_out, VM_out
    
    
    
def v_renegotiated_loop(setup,V,kappa=0.45,interpolate=False):
        # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    
    ism, psm = transition_uniform(setup.agrid,kappa*setup.agrid)
    isf, psf = transition_uniform(setup.agrid,kappa*setup.agrid)
    ind, izf, izm, ipsi = setup.all_indices()
    VMval_single, VFval_single = V['Male, single']['V'], V['Female, single']['V']
    Vval_postren, VMval_postren, VFval_postren = V['Couple']['V'], V['Couple']['VM'], V['Couple']['VF']
    Vm_divorce = (VMval_single[ism,:]*psm[:,None] + VMval_single[ism+1,:]*(1-psm[:,None]))[:,izm]
    Vf_divorce = (VFval_single[isf,:]*psf[:,None] + VFval_single[isf+1,:]*(1-psf[:,None]))[:,izf]
    tg = setup.thetagrid
    
    if interpolate:
        return v_ren_core_int(VFval_postren,VMval_postren,Vval_postren,Vf_divorce,Vm_divorce,tg)
    else:
        return v_ren_core_noint(VFval_postren,VMval_postren,Vval_postren,Vf_divorce,Vm_divorce,tg)

    
    
    
#@jit(nopython=True)#,parallel=True)
def v_ren_core_noint(VF_m,VM_m,V_m,VF_s,VM_s,thetagrid):
    
    na = VF_m.shape[0]
    nz = VF_m.shape[1]
    nt = VF_m.shape[2]
    
    
    VF_out, VM_out, V_out = np.copy(VF_m),  np.copy(VM_m), np.copy(V_m)
    
    for ia in range(na):
        for iz in range(nz):
            vf   = VF_m[ia,iz,:]
            vf_d = VF_s[ia,iz]
            vm   = VM_m[ia,iz,:]
            vm_d = VM_s[ia,iz]
            vc   =  V_m[ia,iz,:]
            
            vf_ren = np.copy(vf)
            vm_ren = np.copy(vm)
            vc_ren = np.copy(vc)
            
            f_agree = ( vf > vf_d )
            m_agree = ( vm > vm_d )
            
            
            for n_f in range(nt):
                    if f_agree[n_f]: break
            
            for n_m in range(nt-1,-1,-1):
                if m_agree[n_m]: break
            
            
            both_agree = (n_f <= n_m)
            
            #both_agree = (f_agree & m_agree)
            
            if not np.any(both_agree):
                # divorce    
                vf_ren[:] = vf_d
                vm_ren[:] = vm_d
                vc_ren[:] = thetagrid*vf_d + (1-thetagrid)*vm_d
            else:
                # renegotiation
                # notice the array slicing rules that are a bit odd
                
                vf_ren[:n_f] = vf_ren[n_f]
                vm_ren[:n_f] = vm_ren[n_f]
                vc_ren[:n_f] = vc_ren[n_f]
                
                vf_ren[n_m+1:] = vf_ren[n_m] 
                vm_ren[n_m+1:] = vm_ren[n_m]
                vc_ren[n_m+1:] = vc_ren[n_m]
            
            VF_out[ia,iz,:] = vf_ren
            VM_out[ia,iz,:] = vm_ren
            V_out [ia,iz,:]  = vc_ren
            
    return V_out, VF_out, VM_out


#@jit(nopython=True)
def v_ren_core_int(VF_m,VM_m,V_m,VF_s,VM_s,thetagrid):
    
    
    na = VF_m.shape[0]
    nz = VF_m.shape[1]
    nt = VF_m.shape[2]
    
    
    NM = np.empty((na,nz),dtype=np.int64)
    NF = np.empty((na,nz),dtype=np.int64)
    
    VF_out, VM_out, V_out = np.copy(VF_m),  np.copy(VM_m), np.copy(V_m)
    
    
    
    for ia in range(na):
        for iz in range(nz):
            vf   = VF_m[ia,iz,:]
            vf_d = VF_s[ia,iz]
            vm   = VM_m[ia,iz,:]
            vm_d = VM_s[ia,iz]
            vc   =  V_m[ia,iz,:]
            
            vf_ren = np.copy(vf)
            vm_ren = np.copy(vm)
            vc_ren = np.copy(vc)
            
            f_agree = ( vf > vf_d )
            m_agree = ( vm > vm_d )
            
            
            for n_f in range(nt):
                if f_agree[n_f]: break
            
            for n_m in range(nt-1,-1,-1):
                if m_agree[n_m]: break
            
            NM[ia,iz] = n_m
            NF[ia,iz] = n_f
            
            
            both_agree  = (n_f <= n_m)
            on_the_edge = (n_f == n_m+1)

            
            if (not both_agree) and (not on_the_edge):
                VF_out[ia,iz,:] = vf_d
                VM_out[ia,iz,:] = vm_d
                V_out [ia,iz,:] = thetagrid*vf_d + (1-thetagrid)*vm_d                
                continue
            
            # here it is possible to reach an agreement
            
            fix_nm, fix_nf = False, False
            
            # this things should not happen but just in case
            if n_m == nt-1: fix_nm, n_m = True, nt-2 # if even at max theta male agrees
            if n_f == 0: fix_nf, n_f = True, 1 # if even at min female agrees
            
            
            a_l = vf[n_f-1] - vf_d
            a_r = vf[n_f] - vf_d    
            kr_f = -a_l/(a_r - a_l) # weight or the right (larger theta) value
            if fix_nf: kr_f = 0.0
            
            b_l = vm[n_m]  - vm_d
            b_r = vm[n_m+1] - vm_d
            kr_m = b_l/(b_l - b_r)  # weight or the right (larger theta) value
            if fix_nm: kr_m = 1.0
            
            
            
            if not fix_nf: assert a_l <= 0 and a_r >= 0
            assert 0 <= kr_f <= 1
            if not fix_nm: assert b_l >= 0 and b_r <= 0  
            
                
            assert 0 <= kr_m <= 1
            
            
            if on_the_edge and kr_f > kr_m:
                # termineate if on the edge and disagreed
                VF_out[ia,iz,:] = vf_d
                VM_out[ia,iz,:] = vm_d
                V_out [ia,iz,:] = thetagrid*vf_d + (1-thetagrid)*vm_d
                continue
                
            
            # renegotiation with linear interpolation
            # notice the array slicing rules that are a bit odd
            
            
            vf_ren_fbind = (1-kr_f)*vf[n_f-1] + kr_f*vf[n_f]
            if not fix_nf: assert np.abs(vf_d - vf_ren_fbind) < 1e-4
            vm_ren_fbind = (1-kr_f)*vm[n_f-1] + kr_f*vm[n_f]
            vc_ren_fbind = (1-kr_f)*vc[n_f-1] + kr_f*vc[n_f]
            
            # at n_f female is already happy
            vf_ren[:n_f] = vf_ren_fbind
            
            
            
            
            vm_ren[:n_f] = vm_ren_fbind
            vc_ren[:n_f] = vc_ren_fbind
            
            
            
            vf_ren_mbind = (1-kr_m)*vf[n_m] + kr_m*vf[n_m+1]
            vm_ren_mbind = (1-kr_m)*vm[n_m] + kr_m*vm[n_m+1]                
            if not fix_nm: assert np.abs(vm_d - vm_ren_mbind) < 1e-4
            vc_ren_mbind = (1-kr_m)*vc[n_m] + kr_m*vc[n_m+1]
            
            # at n_m male is already happy
            vf_ren[n_m+1:] = vf_ren_mbind
            vm_ren[n_m+1:] = vm_ren_mbind
            vc_ren[n_m+1:] = vc_ren_mbind
            
            assert np.all(vf_ren >= vf_d - 1e-4)
            assert np.all(vm_ren >= vm_d - 1e-4)
            
            VF_out[ia,iz,:] = vf_ren
            VM_out[ia,iz,:] = vm_ren
            V_out [ia,iz,:] = vc_ren
            

    
    return V_out, VF_out, VM_out
