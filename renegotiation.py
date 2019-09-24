#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This collects things relevant for renegotiation of couples
"""
from trans_unif import transition_uniform
import numpy as np
from aux_routines import first_true, last_true


def v_last_period_renegotiated(setup,V,kappa=0.45,return_all=False):
    # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    ism, psm = transition_uniform(setup.agrid,kappa*setup.agrid)
    isf, psf = transition_uniform(setup.agrid,kappa*setup.agrid)
    
    ind, izf, izm, ipsi = setup.all_indices()
    
    
    VMval_single, VFval_single = V['SM']['V'], V['SF']['V']
    VMval_postren, VFval_postren = V['M']['VM'], V['M']['VF']
    
    
    Vm_divorce = (VMval_single[ism,:]*psm[:,None] + VMval_single[ism+1,:]*(1-psm[:,None]))[:,izm,None]
    Vf_divorce = (VFval_single[isf,:]*psf[:,None] + VFval_single[isf+1,:]*(1-psf[:,None]))[:,izf,None]
    
    S_f = VFval_postren - Vf_divorce # surplus of female
    S_m = VMval_postren - Vm_divorce # surplus of male
    
    I_f = np.array(S_f > 0) # whether female agrees at this gridpoint
    I_m = np.array(S_m > 0) # whether male agrees at this gridpoint
    
    sq = (I_f & I_m)
    
    together = np.any(sq,axis=2)
    f_ren = (~I_f & together[:,:,None])
    m_ren = (~I_m & together[:,:,None])
    
    nf = first_true(I_f,axis=2)
    nm = last_true(I_m,axis=2)
    
    
    
    # debugging
    Sf_min  = np.take_along_axis(S_f,nf[:,:,None],2).squeeze()
    Sf_min1 = np.take_along_axis(S_f,nf[:,:,None]-1,2).squeeze()
    Sm_min  = np.take_along_axis(S_m,nm[:,:,None],2).squeeze()
    i_sm = np.minimum(nm[:,:,None]+1, setup.ntheta-1) # o/w can exceed the size
    Sm_min1 = np.take_along_axis(S_m,i_sm,2).squeeze()
    assert np.all(Sf_min[together]>0)
    assert np.all(Sf_min1[together]<0)
    assert np.all(Sm_min[together]>0)
    assert np.all(Sm_min1[together]<0)
    # end debugging
    
    # I create new for preren
    VF_out = np.copy(VFval_postren)
    VM_out = np.copy(VMval_postren)
    
    
    Vf_ren_f = np.take_along_axis(VFval_postren,nf[:,:,None],2)
    Vf_ren_m = np.take_along_axis(VFval_postren,nm[:,:,None],2)
    Vm_ren_f = np.take_along_axis(VMval_postren,nf[:,:,None],2)
    Vm_ren_m = np.take_along_axis(VMval_postren,nm[:,:,None],2)
    
    
    bt = lambda x : np.broadcast_to(x, VF_out.shape) # mad skillz
    # this assumed VF_out and VM_out have the same shape
    
    bool_divorce = bt(~together[:,:,None])
    VF_out[bool_divorce] = bt(Vf_divorce)[bool_divorce]
    VM_out[bool_divorce] = bt(Vm_divorce)[bool_divorce]
    
    VF_out[f_ren] = bt(Vf_ren_f)[f_ren]
    VM_out[f_ren] = bt(Vm_ren_f)[f_ren]
    VF_out[m_ren] = bt(Vf_ren_m)[m_ren]
    VM_out[m_ren] = bt(Vm_ren_m)[m_ren]
    
    if not return_all:
        return VF_out, VM_out
    else:
        return VF_out, VM_out, (S_f, S_m, together, f_ren, m_ren, nf, nm)