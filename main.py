#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

import numpy as np
import sobol_seq
from scipy.stats import norm

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
beta = 0.95

# let's approximate three Markov chains


from timeit import default_timer

start = default_timer()

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


def Vlast(sm,sf,zm,zf,psi,theta):
    income = sm + np.exp(zm) + sf + np.exp(zf)
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

def Vs(s,z):    
    return utility(s+np.exp(z))
    
    
#out = Vlast(0,0,grid[-1][20,0],grid[-1][20,1],grid[-1][20,2],0.5)

gamma = 0.5



print('Time elapsed is {}'.format(default_timer()-start))



def theta_opt(sm,sf,zm,zf,psi,npoints=100):
    
    # everything should be an array
    zm = zm if isinstance(zm,np.ndarray) else np.array([zm])
    zf = zf if isinstance(zf,np.ndarray) else np.array([zf])
    psi = psi if isinstance(psi,np.ndarray) else np.array([psi])
    sm = sm if isinstance(sm,np.ndarray) else np.array([sm])
    sf = sf if isinstance(sf,np.ndarray) else np.array([sf])
    
    VM_S = Vs(sm,zm)
    VF_S = Vs(sf,zf)
    
    thetas = np.linspace(0.001,0.999,npoints)
    
    V, VM, VF = Vlast(sm[:,np.newaxis],sf[:,np.newaxis],zm[:,np.newaxis],zf[:,np.newaxis],psi[:,np.newaxis],thetas)
    S_M = (VM - VM_S[:,np.newaxis])
    S_F = (VF - VF_S[:,np.newaxis])
    
    i_pos = (np.array(S_M > 0) & np.array(S_F > 0))
    ind_no = np.where(~i_pos)
    ind_y  =  np.where(i_pos)
    
    nbs_all = np.empty_like(S_M)
    
    nbs_all[ind_no] = -np.inf
    nbs_all[ind_y] = (S_M[ind_y]**gamma)*(S_F[ind_y]**(1-gamma))
    
    # this is matrix containing Nash Bargaining Surplus where rows are states
    # and columns are thetas
    
    # maximize each row over thetas and pick theta    
    theta = np.array([
                        thetas[np.argmax(nbs_row)] if np.any(nbs_row>0)
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
    theta[ind_y]  = thetas[np.argmax(nbs_all[ind_y,:],axis=1)]
    '''
        
    return theta



def Vnext(sm,sf,zm,zf,psi):
    topt = theta_opt(sm,sf,zm,zf,psi)
    
    if sm.size == 1: sm = sm*np.ones_like(sf)
    if sf.size == 1: sf = sf*np.ones_like(sm)
    if zm.size == 1: zm = zm*np.ones_like(zf)
    if zf.size == 1: zf = zf*np.ones_like(zm)
    
    Vout_f = np.zeros_like(zf)
    Vout_m = np.zeros_like(zm)
    Vf = Vs(sf,zf)
    Vm = Vs(sm,zm)
    
    i_mar = ~np.isnan(topt)
    Vmar_tuple = Vlast(sm[i_mar],sf[i_mar],zm[i_mar],zf[i_mar],psi[i_mar],topt[i_mar])
    Vout_m[i_mar] = Vmar_tuple[1]
    Vout_f[i_mar] = Vmar_tuple[2]
    Vout_m[~i_mar] = Vm[~i_mar]
    Vout_f[~i_mar] = Vf[~i_mar]
    
    return Vout_f, Vout_m
    





npart = 7
v = norm.ppf(sobol_seq.i4_sobol_generate(3,npart))
sig_a = 0.1
sig_z = 0.2


def EVnext_f(sf,zf):
    sf = sf if isinstance(sf,np.ndarray) else np.array([sf])
    zf = zf if isinstance(zf,np.ndarray) else np.array([zf])
    sm = sf*np.exp(sig_a*v[:,0])
    zm = zf + sig_z*v[:,1]
    psi = sig_psi_0*v[:,2]
    EV_f = Vnext(sm,sf,zm,zf,psi)[0].mean()
    return EV_f
    
    
from scipy.optimize import fminbound


nnodes = 7
z_nodes = norm.ppf(sobol_seq.i4_sobol_generate(1,nnodes))

def Vzero(a0,z0):
    income = a0 + np.exp(z0)
    z_next = z0 + sig_zf*z_nodes
    
    def total_u(s):
        c = income - s
        
        EV = np.mean([EVnext_f(s,z) for z in z_next])
        return -(utility(c) + beta*EV)
    
    
    smin = 0
    smax = 0.9*income
    
    s_opt = fminbound(total_u,smin,smax)
    
    return -total_u(s_opt), s_opt, s_opt/income


#print(Vzero(0,0.4))
        
    
    

    
    
#tht = theta_opt(0,0,grid[-1][:,0],grid[-1][:,1],grid[-1][:,2])

#Vnxt = Vnext(0,0,grid[-1][:,0],grid[-1][:,1],grid[-1][:,2])
   

a0v = np.arange(0,5,0.01)

savings = [Vzero(a0,-0.4)[2] for a0 in a0v]

import matplotlib.pyplot as plt
plt.figure()
plt.plot(a0v,savings)


#print(tht)
print('Time elapsed is {}'.format(default_timer()-start))
