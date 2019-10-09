#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects routines for renegotiation and marriage

Some code is meant to be reused

"""


from trans_unif import transition_uniform
import numpy as np
from aux_routines import first_true, last_true


# this is renegotiator
def v_ren(setup,V,kappa=0.45,interpolate=False):
    # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    ism, psm = transition_uniform(setup.agrid,kappa*setup.agrid)
    isf, psf = transition_uniform(setup.agrid,kappa*setup.agrid)
    
    ind, izf, izm, ipsi = setup.all_indices()
    
    
    #nt = setup.thetagrid.size
    
    VMval_single, VFval_single = V['Male, single']['V'], V['Female, single']['V']
    Vval_postren, VMval_postren, VFval_postren = V['Couple']['V'], V['Couple']['VM'], V['Couple']['VF']
    
    
    Vm_divorce = (VMval_single[ism,:]*psm[:,None] + VMval_single[ism+1,:]*(1-psm[:,None]))[:,izm]
    Vf_divorce = (VFval_single[isf,:]*psf[:,None] + VFval_single[isf+1,:]*(1-psf[:,None]))[:,izf]
    
    
    return v_ren_core(VFval_postren,VMval_postren,Vval_postren,Vf_divorce,Vm_divorce,setup.thetagrid,interpolate=interpolate)


def v_ren_core(VF_yes,VM_yes,VC_yes,VF_no,VM_no,thetagrid,interpolate=False):
    
    
    
    ax = VC_yes.ndim - 1 # axis where theta is
    nt = VC_yes.shape[ax] # number of thetas
    
    assert nt==thetagrid.size
    
    print((ax,nt))
    
    shp = VF_yes.shape[:-1] + (1,)
    
    VF_no = VF_no.reshape(shp)
    VM_no = VM_no.reshape(shp)
    
    S_f = VF_yes - VF_no # surplus of female
    S_m = VM_yes - VM_no # surplus of male
    
    I_f = np.array(S_f > 0) # whether female agrees at this gridpoint
    I_m = np.array(S_m > 0) # whether male agrees at this gridpoint
    
    
    sq = (I_f & I_m)
    
    
    nf = first_true(I_f,axis=ax)
    nm = last_true(I_m,axis=ax)
    
    together = np.any(sq,axis=ax)
    
    
    
    
    together_2 = (np.array(nf<=nm) & np.array(nf!=-1) & np.array(nm!=-1) )
    
    on_the_edge = (np.array(nf == nm+1) & np.array(nf!=-1) & np.array(nm!=-1) ).reshape(shp)
    
    assert np.all(together_2 == together)
    
    
    NF = nf.reshape(shp)
    NF_fix = np.array( NF == 0 ) # when first true is zero
    NF[NF_fix] = 1
    NF_ok = np.array( NF != -1 ) # when there is at least some true value
    NF_sc = (NF_ok & ~NF_fix) # whether it corresponds to sign change
    NF[~NF_ok] = nt-1
    
    
    
    NM = nm.reshape(shp)
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
    
    VF_rf, VF_lf, VF_rm, VF_lm  = (np.take_along_axis(VF_yes,x,ax) for x in (rf, lf, rm, lm))
    VM_rf, VM_lf, VM_rm, VM_lm =  (np.take_along_axis(VM_yes,x,ax) for x in (rf, lf, rm, lm))
    VC_rf, VC_lf, VC_rm, VC_lm =  (np.take_along_axis( VC_yes,x,ax) for x in (rf, lf, rm, lm))
    

    A_l, A_r = VF_lf - VF_no, VF_rf - VF_no
    B_l, B_r = VM_lm - VM_no, VM_rm - VM_no
    
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
    
    
    
    
    
    
    # the further part is specific to renegotiation
    

 
    
    f_ren = (~I_f & together.reshape(shp)) # female-initiated renegotiation
    m_ren = (~I_m & together.reshape(shp)) # male-initiated renegotiation
     
    
    
    # I create new for preren
    VF_out = np.copy(VF_yes)
    VM_out = np.copy(VM_yes)
    V_out  = np.copy(VC_yes)
    
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
        
        
        
        
        
        assert np.all( np.abs( Vf_ifren_f - VF_no)[NF_sc] < 1e-4 )
        assert np.all( np.abs( Vm_ifren_m - VM_no)[NM_sc] < 1e-4 )
        
        
        
    bt = lambda x : np.broadcast_to(x, VF_out.shape) # mad skillz
    # this assumed VF_out and VM_out have the same shape
    
    tshape = (VF_yes.ndim-1)*(1,) + (nt,)
    
    t_stretch = bt(thetagrid.reshape(tshape))
    
    
    
    bool_divorce = bt(~together.reshape(shp))
    VF_out[bool_divorce] = bt(VF_no)[bool_divorce]
    VM_out[bool_divorce] = bt(VM_no)[bool_divorce]
    
    V_out[bool_divorce] = t_stretch[bool_divorce]*VF_out[bool_divorce] + \
                      (1-t_stretch[bool_divorce])*VM_out[bool_divorce] # check me please
    
    VF_out[f_ren] = bt(Vf_ifren_f)[f_ren]
    VM_out[f_ren] = bt(Vm_ifren_f)[f_ren]
    V_out[f_ren]  = bt( V_ifren_f)[f_ren]
    
    VF_out[m_ren] = bt(Vf_ifren_m)[m_ren]
    VM_out[m_ren] = bt(Vm_ifren_m)[m_ren]
    V_out[m_ren]  = bt( V_ifren_m)[m_ren]
    
    
    
    
    assert np.all(VF_out >= VF_no - 1e-4)
    assert np.all(VM_out >= VM_no - 1e-4)
    
    return V_out, VF_out, VM_out