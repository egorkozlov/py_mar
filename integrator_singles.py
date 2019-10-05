#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines for intergation for singles
"""

import numpy as np
from trans_unif import transition_uniform
from mc_tools import int_prob


def v_after_mar_grid(setup,V,sf,sm,ind_or_inds):
    # this is generinc and does not depend on gender
    
    # setting up
    agrid = setup.agrid
    gamma = setup.pars['m_bargaining_weight']    
    VMval_single, VFval_single = V['Male, single']['V'], V['Female, single']['V']
    VMval_postren, VFval_postren = V['Couple']['VM'], V['Couple']['VF']
    
    
    # substantial part
    ind, izf, izm, ipsi = setup.all_indices(ind_or_inds)
    
    isf, psf = transition_uniform(agrid,sf)
    ism, psm = transition_uniform(agrid,sm)
    isc, psc = transition_uniform(agrid,sf+sm)
    
    #assert np.all(setup.exogrid.all_t[-1][ind,0] == setup.exogrid.zf_t[-1][izf])
    #assert np.all(setup.exogrid.all_t[-1][ind,1] == setup.exogrid.zm_t[-1][izm])
    #assert np.all(setup.exogrid.all_t[-1][ind,2] == setup.exogrid.psi_t[-1][ipsi])
    
    ism, isf, isc, psf, psm, psc = (x[:,None] for x in (ism,isf,isc,psf, psm, psc))
    
    
    Vms = VMval_single[ism,izm]*psm + VMval_single[ism+1,izm]*(1-psm)
    Vfs = VFval_single[isf,izf]*psf + VFval_single[isf+1,izf]*(1-psf)
    Vmm = VMval_postren[isc,ind,:]*(psc[:,:,None]) + VMval_postren[isc+1,ind,:]*(1-psc[:,:,None])
    Vfm = VFval_postren[isc,ind,:]*(psc[:,:,None]) + VFval_postren[isc+1,ind,:]*(1-psc[:,:,None])
    
    s_m = Vmm - Vms[:,:,None]
    s_f = Vfm - Vfs[:,:,None]
    
    nbs = -np.inf*np.ones_like(s_m)
    
    i_pos = ((s_m > 0) & (s_f>0) )
    nbs[i_pos] = s_m[i_pos]**(gamma) * s_f[i_pos]**(1-gamma)
    
    ismar = np.any(nbs>0,axis=2)
    Vout_m, Vout_f = np.empty_like(Vms), np.empty_like(Vfs)
    
    Vout_m[~ismar] = Vms[~ismar]
    Vout_m[ismar]  = Vmm[ismar,nbs[ismar,:].argmax(axis=1)]
    Vout_f[~ismar] = Vfs[~ismar]
    Vout_f[ismar]  =  Vfm[ismar,nbs[ismar,:].argmax(axis=1)]
    
    
    return Vout_f, Vout_m
    

def ev_after_savings_grid_all_z(setup,V,sown,female,t,trim_lvl=0.01):
    # this takes gender as argument so should be called twice
    
    nexo = setup.pars['nexo']
    sigma_psi_init = setup.pars['sigma_psi_init']
    sig_z_partner = setup.pars['sig_partner_z']
    sig_a_partner = setup.pars['sig_partner_a']
    ns = sown.size
    eps_a_partner = setup.integration['nodes_couple'][:,0]
    npart = setup.integration['num_partners']    
    psi_couple = setup.exogrid.psi_t[-1]
    
    
    
    if female:
        nz_single = setup.exogrid.zf_t[t].shape[0]
        p_mat = np.empty((nexo,nz_single))
        z_own = setup.exogrid.zf_t[t]
        n_zown = z_own.shape[0]
        z_partner = setup.exogrid.zm_t[t+1]
        zmat_own = setup.exogrid.zf_t_mat[t]
        i_vnext = 0
    else:
        nz_single = setup.exogrid.zm_t[t].shape[0]
        p_mat = np.empty((nexo,nz_single))
        z_own = setup.exogrid.zm_t[t]
        n_zown = z_own.shape[0]
        z_partner = setup.exogrid.zf_t[t+1]
        zmat_own = setup.exogrid.zm_t_mat[t]    
        i_vnext = 1
        
    
    def ind_conv(a,b,c): return setup.all_indices((a,b,c))[0]
    
    
    for iz in range(n_zown):
        p_psi = int_prob(psi_couple,mu=0,sig=sigma_psi_init)
        if female:
            p_zm  = int_prob(z_partner, mu=z_own[iz],sig=sig_z_partner)
            p_zf  = zmat_own[iz,:]
        else:
            p_zf  = int_prob(z_partner, mu=z_own[iz],sig=sig_z_partner)
            p_zm  = zmat_own[iz,:]
        #sm = sf
    
        p_vec = np.zeros(nexo)
        
        for izf, p_zf_i in enumerate(p_zf):
            for izm, p_zm_i in enumerate(p_zm):
                for ipsi, p_psi_i in enumerate(p_psi):
                    p = p_zf_i*p_zm_i*p_psi_i
                    if p > trim_lvl:
                        p_vec[ind_conv(izf,izm,ipsi)] = p    
                        
                        
        p_vec = p_vec / np.sum(p_vec)
        p_mat[:,iz] = p_vec
        
    
    V_next = np.ones((ns,nexo))*(-1e-10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    EV = 0.0
    
    for ipart in range(npart):
        sm = sown*np.exp(sig_a_partner*eps_a_partner[ipart])
        V_next[:,inds] = v_after_mar_grid(setup,V,sown,sm,inds)[i_vnext]
        EV += (1/npart)*np.dot(V_next,p_mat)
        
    
    return EV#, p_mat, VF_next