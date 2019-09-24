#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fminbound
from timeit import default_timer
from setup import ModelSetup
from aux_routines import first_true, last_true


plt.figure()
start = default_timer()


setup = ModelSetup()



print('Time elapsed is {}'.format(default_timer()-start))


# grid for s

gexo = setup.exogrid.all_t[-1]
nexo = gexo.shape[0]


# relevant for plotting
a0 = 2
zgrid = setup.exogrid.zf_t[0]
z0 = setup.exogrid.zf_t[0][4]




  
#out = setup.vm_last(0,0,setup.exogrid[-1][20,0],setup.exogrid[-1][20,1],setup.exogrid[-1][20,2],0.5)



print('Time elapsed is {}'.format(default_timer()-start))
#raise Exception('stop here')


# two auxiliary routines to parse different kinds of inputs


from main_exact import exact_solution, exact_solution_with_interpolated_ev


savings_rate = exact_solution(setup,setup.agrid,np.array([z0]),ireturn=2)

plt.plot(setup.agrid,savings_rate,label="exact") # these are exact savings on the grid


#print(tht)
print('Time elapsed is {}'.format(default_timer()-start))


savings_rate_int = exact_solution_with_interpolated_ev(setup,setup.agrid,np.array([z0]),ireturn=2)
plt.plot(setup.agrid,savings_rate_int,label="interpolation-dumb")


print('Time elapsed is {}'.format(default_timer()-start))


## this section uses grids for everything
def vm_last_grid():
    Vval_postren, VMval_postren, VFval_postren =  \
    np.zeros((setup.na,setup.nexo,setup.ntheta)), np.zeros((setup.na,setup.nexo,setup.ntheta)), np.zeros((setup.na,setup.nexo,setup.ntheta))

    for itheta in range(setup.ntheta):
        Vval_postren[:,:,itheta], VMval_postren[:,:,itheta], VFval_postren[:,:,itheta] = setup.vm_last(setup.agrid[:,None],gexo[:,0],gexo[:,1],gexo[:,2],setup.thetagrid[itheta])
    return Vval_postren, VMval_postren, VFval_postren

Vval_postren, VMval_postren, VFval_postren = vm_last_grid()

VMval_single, VFval_single = setup.vs_last(setup.agrid[:,None],setup.exogrid.zm_t[-1]), setup.vs_last(setup.agrid[:,None],setup.exogrid.zf_t[-1])

print('Time elapsed is {}'.format(default_timer()-start))


from trans_unif import transition_uniform
from mc_tools import int_prob



def v_after_mar_grid(sf,sm,ind_or_inds):
    
    ind, izf, izm, ipsi = setup.all_indices(ind_or_inds)
    
    isf, psf = transition_uniform(setup.agrid,sf)
    ism, psm = transition_uniform(setup.agrid,sm)
    isc, psc = transition_uniform(setup.agrid,sf+sm)
    
    assert np.all(gexo[ind,0] == setup.exogrid.zf_t[-1][izf])
    assert np.all(gexo[ind,1] == setup.exogrid.zm_t[-1][izm])
    assert np.all(gexo[ind,2] == setup.exogrid.psi_t[-1][ipsi])
    
    ism, isf, isc, psf, psm, psc = (x[:,None] for x in (ism,isf,isc,psf, psm, psc))
    
    
    Vms = VMval_single[ism,izm]*psm + VMval_single[ism+1,izm]*(1-psm)
    Vfs = VFval_single[isf,izf]*psf + VFval_single[isf+1,izf]*(1-psf)
    Vmm = VMval_postren[isc,ind,:]*(psc[:,:,None]) + VMval_postren[isc+1,ind,:]*(1-psc[:,:,None])
    Vfm = VFval_postren[isc,ind,:]*(psc[:,:,None]) + VFval_postren[isc+1,ind,:]*(1-psc[:,:,None])
    
    s_m = Vmm - Vms[:,:,None]
    s_f = Vfm - Vfs[:,:,None]
    
    nbs = -np.inf*np.ones_like(s_m)
    
    i_pos = ((s_m > 0) & (s_f>0) )
    nbs[i_pos] = s_m[i_pos]**(setup.pars['m_bargaining_weight']) * s_f[i_pos]**(1-setup.pars['m_bargaining_weight'])
    
    ismar = np.any(nbs>0,axis=2)
    Vout_m, Vout_f = np.empty_like(Vms), np.empty_like(Vfs)
    
    Vout_m[~ismar] = Vms[~ismar]
    Vout_m[ismar]  = Vmm[ismar,nbs[ismar,:].argmax(axis=1)]
    Vout_f[~ismar] = Vfs[~ismar]
    Vout_f[ismar]  =  Vfm[ismar,nbs[ismar,:].argmax(axis=1)]
    
    
    return Vout_f, Vout_m#, Vms, Vfs, Vmm, Vfm, nbs
    

def ev_after_savings_grid_all_z(sf,trim_lvl=0.01):
    
    p_mat = np.empty((setup.pars['nexo'],setup.exogrid.zf_t[0].shape[0]))
    
    
    for iz in range(setup.exogrid.zf_t[0].shape[0]):
        p_zm  = int_prob(setup.exogrid.zm_t[-1], mu=setup.exogrid.zf_t[0][iz],sig=setup.pars['sig_partner_z'])
        p_psi = int_prob(setup.exogrid.psi_t[-1],mu=0,sig=setup.pars['sigma_psi_init'])
        p_zf  = setup.exogrid.zf_t_mat[0][iz,:]
        #sm = sf
    
        p_vec = np.zeros(setup.pars['nexo'])
        
        
        
        for izf, p_zf_i in enumerate(p_zf):
            for izm, p_zm_i in enumerate(p_zm):
                for ipsi, p_psi_i in enumerate(p_psi):
                    p = p_zf_i*p_zm_i*p_psi_i
                    if p > trim_lvl:
                        p_vec[izf*setup.pars['n_zm']*setup.pars['n_psi'] + izm*setup.pars['n_psi'] + ipsi] = p
        
        p_vec = p_vec / np.sum(p_vec)
        p_mat[:,iz] = p_vec
        
    
    VF_next = np.ones((sf.size,setup.pars['nexo']))*(-1e-10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    EV_f = 0.0
    
    for ipart in range(setup.integration['num_partners']):
        sm = sf*np.exp(setup.pars['sig_partner_a']*setup.integration['nodes_couple'][ipart,0])
        VF_next[:,inds] = v_after_mar_grid(sf,sm,inds)[0]
        EV_f += (1/setup.integration['num_partners'])*np.dot(VF_next,p_mat)
        
    
    
    
    return EV_f#, p_mat, VF_next


EV_integrated = ev_after_savings_grid_all_z(setup.agrid)


def v_period_zero_grid(a0):
    income = a0[:,None] + np.exp(setup.exogrid.zf_t[0][None,:])
    
    
    def neg_total_u(s,inc,EVg):
        c = inc - s
        assert np.all(c > 0)
        EV = np.interp(s,setup.agrid,EVg)        
        return -(setup.u(c) + setup.pars['beta']*EV)
    
    
    
    smin = np.zeros_like(income)
    smax = 0.9*income
    
    EVT = EV_integrated.T
    
    s_opt = np.array([
                         [
                          fminbound( lambda x : neg_total_u(x,income[j,i],EVval),smin[j,i],smax[j,i] ) 
                          for i, EVval in enumerate(EVT)
                         ]
                       for j in range(a0.size)
                      ])
    
    
    V_ret =  np.array([-neg_total_u(s_opt[:,i],income[:,i],EV) for i, EV in enumerate(EVT)]).T
    
    return V_ret, s_opt, s_opt/income



def v_last_period_renegotiated(kappa=0.45,return_all=False):
    # this returns value functions for couple that entered the last period with
    # (s,Z,theta) from the grid and is allowed to renegotiate them or breakup
    
    ism, psm = transition_uniform(setup.agrid,kappa*setup.agrid)
    isf, psf = transition_uniform(setup.agrid,kappa*setup.agrid)
    
    ind, izf, izm, ipsi = setup.all_indices()
    
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

vv = v_last_period_renegotiated()
    

print('Time elapsed is {}'.format(default_timer()-start))
sss = v_period_zero_grid(setup.agrid)
plt.plot(setup.agrid,sss[2][:,4],label="interpolation-smart")
plt.legend()
print('Time elapsed is {}'.format(default_timer()-start))


