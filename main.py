#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

import numpy as np

T         = 2
sig_zf_0  = 0.15
sig_zf    = 0.05
n_zf      = 9
sig_zm_0  = 0.2
sig_zm    = 0.075
n_zm      = 5
sig_psi_0 = 0.12
sig_psi   = 0.03
n_psi     = 15

# let's approximate three Markov chains


from rw_approximations import rouw_nonst
from mc_tools import combine_matrices_two_lists


zf,   zf_mat = rouw_nonst(T,sig_zf,sig_zf_0,n_zf)
zm,   zm_mat = rouw_nonst(T,sig_zm,sig_zm_0,n_zm)
psi, psi_mat = rouw_nonst(T,sig_psi,sig_psi_0,n_psi)


zfzm, zfzm_mat = combine_matrices_two_lists(zf,zm,zf_mat,zm_mat)
grid, mat = combine_matrices_two_lists(zfzm,psi,zfzm_mat,psi_mat)

A = 1.2
sigma = 1.5
rho = 0.0

def umult(theta):
    assert np.all(theta > 0) and np.all(theta < 1)
    powr = (1+rho)/(rho+sigma)
    tf = theta
    tm = 1-theta
    ces = (tf**powr + tm**powr)**(1/powr)
    return (A**(1-sigma))*ces

#print(umult(0.5))

def cmult(theta):
    assert np.all(theta > 0) and np.all(theta < 1)
    powr = (1+rho)/(rho+sigma)
    irho = 1/(1+rho)
    irs  = 1/(rho+sigma)
    tf = theta
    tm = 1-theta
    bottom = (tf**(powr) + tm**(powr))**irho  
    
    kf = A*(tf**(irs))/bottom
    km = A*(tm**(irs))/bottom
    return kf, km

def utility(c):
    return (c**(1-sigma))/(1-sigma)

#print(cmult(0.5))


def Vlast(zm,zf,psi,theta):
    income = np.exp(zm) + np.exp(zf)
    kf, km = cmult(theta)
    cf, cm = kf*income, km*income
    u_couple = umult(theta)*utility(income)
    #u_couple_2 = theta*utility(cf) + (1-theta)*utility(cm)
    #print((u_couple-u_couple_2))
    u_m = utility(cm)
    u_f = utility(cf)
    V = u_couple + psi
    VM = u_m + psi
    VF = u_f + psi
    return V, VM, VF

def Vs(z):
    return utility(np.exp(z))
    
    
out = Vlast(grid[-1][20,0],grid[-1][20,1],grid[-1][20,2],0.5)

gamma = 0.5

def theta_opt(zm,zf,psi):
    
    VM_S = Vs(zm)
    VF_S = Vs(zf)
    
    def nbs(theta,ret_s=False):
        V, VM, VF = Vlast(zm,zf,psi,theta)
        S_M = VM - VM_S
        S_F = VF - VF_S
        
        out = -np.inf
        if S_M > 0 and S_F > 0:
            out = (S_M**gamma)*(S_F**(1-gamma))
            
        return out
    
    thetas = np.linspace(0.05,0.95,50)
    
    nbs_all = [nbs(t) for t in thetas]
    
    if np.any(np.array(nbs_all)>0):
        theta = thetas[np.argmax(nbs_all)]
    else:
        theta = np.nan
        
    return theta

tht = np.array([theta_opt(g[0],g[1],g[2]) for g in grid[-1]])



      
print(tht)
