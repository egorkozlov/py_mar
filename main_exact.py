#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This contains routines relevant for obtaining exact solution
"""

import numpy as np
from aux_routines import num_to_nparray, unify_sizes
from scipy.optimize import fminbound


def theta_opt(setup,sm,sf,zm,zf,psi):
    # this finds optimal theta using discrete choice over grid of theta
    
    zm,zf,psi,sm,sf = num_to_nparray(zm,zf,psi,sm,sf)
    
    
    VM_S = setup.vs_last(sm,zm)
    VF_S = setup.vs_last(sf,zf)
    
    
    V, VM, VF = setup.vm_last(sm[:,np.newaxis]+sf[:,np.newaxis],zm[:,np.newaxis],zf[:,np.newaxis],psi[:,np.newaxis],setup.thetagrid)
    S_M = (VM - VM_S[:,np.newaxis])
    S_F = (VF - VF_S[:,np.newaxis])
    
    i_pos = (np.array(S_M > 0) & np.array(S_F > 0))
    ind_no = np.where(~i_pos)
    ind_y  =  np.where(i_pos)
    
    nbs_all = np.empty_like(S_M)
    
    nbs_all[ind_no] = -np.inf
    nbs_all[ind_y] = (S_M[ind_y]**setup.pars['m_bargaining_weight'])*(S_F[ind_y]**(1-setup.pars['m_bargaining_weight']))
    
    # this is matrix containing Nash Bargaining Surplus where rows are states
    # and columns are thetas
    
    # maximize each row over thetas and pick theta    
    theta = np.array([
                        setup.thetagrid[np.argmax(nbs_row)] if np.any(nbs_row>0)
                        else np.nan 
                        for nbs_row in nbs_all]
                    )
      

    ''' 
    # alternative version w/o looping:
    
    i_no = np.all(nbs_all<0,axis=1)
    ind_no = np.where(i_no)[0]
    ind_y  = np.where(~i_no)[0]
    theta = np.empty_like(zm)
    
    theta[ind_no] = np.nan
    theta[ind_y]  = setup.thetagrid[np.argmax(nbs_all[ind_y,:],axis=1)]
    '''
        
    return theta




def v_after_mar_exact(setup,sm,sf,zm,zf,psi):
    
    sm, sf, zm, zf, psi = num_to_nparray(sm,sf,zm,zf,psi)
    
    topt = theta_opt(setup,sm,sf,zm,zf,psi)  
    
    
    sm, sf = unify_sizes(sm,sf)
    zm, zf = unify_sizes(zm,zf)
    
    
    Vout_f = np.zeros_like(zf)
    Vout_m = np.zeros_like(zm)
    Vf = setup.vs_last(sf,zf)
    Vm = setup.vs_last(sm,zm)
    
    i_mar = ~np.isnan(topt)
    
    if np.any(i_mar):
        Vmar_tuple = setup.vm_last(sm[i_mar]+sf[i_mar],zm[i_mar],zf[i_mar],psi[i_mar],topt[i_mar])
        Vout_m[i_mar] = Vmar_tuple[1]
        Vout_f[i_mar] = Vmar_tuple[2]
    
    Vout_m[~i_mar] = Vm[~i_mar]
    Vout_f[~i_mar] = Vf[~i_mar]
    
    return Vout_f, Vout_m
    




def ev_after_savings_exact(setup,sf,zf):
    # only one thing can be vectorized
    # this always returns a one-dimensional vector
    sf = sf if isinstance(sf,np.ndarray) else np.array([sf])
    zf = zf if isinstance(zf,np.ndarray) else np.array([zf])
    
    v = setup.integration['nodes_couple']
    
    V_f = np.zeros((sf.size,v.shape[0]))
    
    for ipart in range(v.shape[0]):
        sm  = sf*np.exp(setup.pars['sig_partner_a']*setup.integration['nodes_couple'][ipart,0])
        zm  = np.ones_like(sf)*(zf + setup.pars['sig_partner_z']*setup.integration['nodes_couple'][ipart,1])
        psi = np.ones_like(sf)*(setup.pars['sigma_psi_init']*setup.integration['nodes_couple'][ipart,2])
        V_f[:,ipart] = v_after_mar_exact(setup,sm,sf,zm,zf,psi)[0]
        
    
    EV_f = V_f.mean(axis=1)
    return EV_f
    
    


def v_period_zero_exact(setup,a0,z0):
    income = a0 + np.exp(z0)
    z_next = z0 + setup.pars['sig_zf']*setup.integration['z_nodes']
    
    def neg_total_u(s):
        c = income - s
        
        EV = np.mean([ev_after_savings_exact(setup,s,z) for z in z_next])
        return -(setup.u(c) + setup.pars['beta']*EV)
    
    smin = 0
    smax = 0.9*income
    
    s_opt = fminbound(neg_total_u,smin,smax)
    
    return -neg_total_u(s_opt), s_opt, s_opt/income



def exact_solution(setup,a_values,z_values,ireturn=1):
    
    # this obtains exact solution for given set of values for a and z
    result =          [
                        [
                           v_period_zero_exact(setup,a,z)[ireturn] 
                            for z in z_values
                        ]
                        for a in a_values
                      ]


    return np.array(result)
    

# this part stretches future value function on grid instead of direct implementation
# this part is pretty useless: it replaces optimization with precomputation of EVnext]


def v_period_zero_interp_at_optimization(setup,a0,z0):
    income = a0 + np.exp(z0)
    z_next = z0 + setup.pars['sig_zf']*setup.integration['z_nodes']
    
    
    EVgrid = np.mean(np.transpose([ev_after_savings_exact(setup,setup.agrid,z) for z in z_next]),axis=1)
    
    
    def neg_total_u(s):
        c = income - s
        EV = np.interp(s,setup.agrid,EVgrid)        
        return -(setup.u(c) + setup.pars['beta']*EV)
    
    
    smin = 0
    smax = 0.9*income
    
    s_opt = fminbound(neg_total_u,smin,smax)
    
    return -neg_total_u(s_opt), s_opt, s_opt/income
    
    
def exact_solution_with_interpolated_ev(setup,a_values,z_values,ireturn=1):    
    # this obtains exact solution for given set of values for a and z
    result =          [
                        [
                           v_period_zero_interp_at_optimization(setup,a,z)[ireturn] 
                            for z in z_values
                        ]
                        for a in a_values
                      ]


    return np.array(result)
    

    