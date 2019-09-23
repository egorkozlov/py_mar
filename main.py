#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""

import numpy as np
import sobol_seq
from scipy.stats import norm
from numba import jit

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

print('Time elapsed is {}'.format(default_timer()-start))

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

#@jit(nopython=True)
def utility(c):
    return (c**(1-sigma))/(1-sigma)

#print(cmult(0.5))

def V_postren(s,zm,zf,psi,theta):
    income = s + np.exp(zm) +  np.exp(zf)
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





ntheta=200
thetas = np.linspace(0.01,0.99,ntheta)
ns = 100
sgrid = np.linspace(0,10,ns)

gexo = grid[-1]

nexo = gexo.shape[0]


#@jit


def Vs(s,z):    
    return utility(s+np.exp(z))
  
#out = V_postren(0,0,grid[-1][20,0],grid[-1][20,1],grid[-1][20,2],0.5)

gamma = 0.5



print('Time elapsed is {}'.format(default_timer()-start))
#raise Exception('stop here')

def theta_opt(sm,sf,zm,zf,psi):
    
    # everything should be an array
    zm = zm if isinstance(zm,np.ndarray) else np.array([zm])
    zf = zf if isinstance(zf,np.ndarray) else np.array([zf])
    psi = psi if isinstance(psi,np.ndarray) else np.array([psi])
    sm = sm if isinstance(sm,np.ndarray) else np.array([sm])
    sf = sf if isinstance(sf,np.ndarray) else np.array([sf])
    
    VM_S = Vs(sm,zm)
    VF_S = Vs(sf,zf)
    
    
    V, VM, VF = V_postren(sm[:,np.newaxis]+sf[:,np.newaxis],zm[:,np.newaxis],zf[:,np.newaxis],psi[:,np.newaxis],thetas)
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
    Vmar_tuple = V_postren(sm[i_mar]+sf[i_mar],zm[i_mar],zf[i_mar],psi[i_mar],topt[i_mar])
    Vout_m[i_mar] = Vmar_tuple[1]
    Vout_f[i_mar] = Vmar_tuple[2]
    Vout_m[~i_mar] = Vm[~i_mar]
    Vout_f[~i_mar] = Vf[~i_mar]
    
    return Vout_f, Vout_m
    



npart = 10
v = norm.ppf(sobol_seq.i4_sobol_generate(3,npart))
sig_a = 0.1
sig_z = 0.2


def EVnext_f(sf,zf):
    # only one thing can be vectorized
    # this always returns a one-dimensional vector
    sf = sf if isinstance(sf,np.ndarray) else np.array([sf])
    zf = zf if isinstance(zf,np.ndarray) else np.array([zf])
    
    V_f = np.zeros((sf.size,v.shape[0]))
    
    for ipart in range(v.shape[0]):
        sm  = sf*np.exp(sig_a*v[ipart,0])
        zm  = np.ones_like(sf)*(zf + sig_z*v[ipart,1])
        psi = np.ones_like(sf)*(sig_psi_0*v[ipart,2])
        V_f[:,ipart] = Vnext(sm,sf,zm,zf,psi)[0]
        
    
    EV_f = V_f.mean(axis=1)
    return EV_f
    
    
from scipy.optimize import fminbound


nnodes = 7
z_nodes = norm.ppf(sobol_seq.i4_sobol_generate(1,nnodes))

def Vzero(a0,z0):
    income = a0 + np.exp(z0)
    z_next = z0 + sig_zf*z_nodes
    
    def neg_total_u(s):
        c = income - s
        
        EV = np.mean([EVnext_f(s,z) for z in z_next])
        return -(utility(c) + beta*EV)
    
    
    smin = 0
    smax = 0.9*income
    
    s_opt = fminbound(neg_total_u,smin,smax)
    
    return -neg_total_u(s_opt), s_opt, s_opt/income


a0 = 2

zgrid = zf[0]


savings_rate = [Vzero(a0,z)[2] for z in zgrid]

import matplotlib.pyplot as plt
plt.figure()
plt.plot(zgrid,savings_rate,label="exact") # these are exact savings on the grid


#print(tht)
print('Time elapsed is {}'.format(default_timer()-start))





# this part stretches future value function on grid instead of direct implementation


# this part is pretty useless: it replaces optimization with precomputation of EVnext
def Vzero_int(a0,z0):
    income = a0 + np.exp(z0)
    z_next = z0 + sig_zf*z_nodes
    
    EVgrid = np.mean(np.transpose([EVnext_f(sgrid,z) for z in z_next]),axis=1)
    
    
    
    def neg_total_u(s):
        c = income - s
        EV = np.interp(s,sgrid,EVgrid)        
        return -(utility(c) + beta*EV)
    
    
    smin = 0
    smax = 0.9*income
    
    s_opt = fminbound(neg_total_u,smin,smax)
    
    return -neg_total_u(s_opt), s_opt, s_opt/income

savings_rate_int = [Vzero_int(a0,z)[2] for z in zgrid]
plt.plot(zgrid,savings_rate_int,label="interpolation")
plt.legend()

print('Time elapsed is {}'.format(default_timer()-start))







## this section uses grids for everything

def Vlast():
    Vval_postren, VMval_postren, VFval_postren =  \
    np.zeros((ns,nexo,ntheta)), np.zeros((ns,nexo,ntheta)), np.zeros((ns,nexo,ntheta))

    for itheta in range(ntheta):
        Vval_postren[:,:,itheta], VMval_postren[:,:,itheta], VFval_postren[:,:,itheta] = V_postren(sgrid[:,None],gexo[:,0],gexo[:,1],gexo[:,2],thetas[itheta])
    return Vval_postren, VMval_postren, VFval_postren

Vval_postren, VMval_postren, VFval_postren = Vlast()

VMval_single, VFval_single = Vs(sgrid[:,None],zm[-1]), Vs(sgrid[:,None],zf[-1])

print('Time elapsed is {}'.format(default_timer()-start))


from trans_unif import transition_uniform
from mc_tools import int_prob


def V_newmar(sf,sm,ind_or_inds):
    
    if isinstance(ind_or_inds,tuple):
        izf,izm,ipsi = ind_or_inds
        ind = izf*n_zm*n_psi + izm*n_psi + ipsi
        #print(ind,izf,izm,ipsi)
    else:
        ind = ind_or_inds
        izf = ind // (n_zm*n_psi)
        izm = (ind - izf*n_zm*n_psi) // n_psi
        ipsi = ind - izf*n_zm*n_psi - izm*n_psi
        #print(ind,izf,izm,ipsi)
    
    
    isf, psf = transition_uniform(sgrid,sf)
    ism, psm = transition_uniform(sgrid,sm)
    isc, psc = transition_uniform(sgrid,sf+sm)
    
    assert np.all(gexo[ind,0] == zf[-1][izf])
    assert np.all(gexo[ind,1] == zm[-1][izm])
    assert np.all(gexo[ind,2] == psi[-1][ipsi])
    
    
    
    Vms = VMval_single[ism,izm]*psm + VMval_single[ism+1,izm]*(1-psm)
    Vfs = VFval_single[isf,izf]*psf + VFval_single[isf+1,izf]*(1-psf)
    Vmm = VMval_postren[isc,ind,:]*psc[:,None] + VMval_postren[isc+1,ind,:]*(1-psc[:,None])
    Vfm = VFval_postren[isc,ind,:]*psc[:,None] + VFval_postren[isc+1,ind,:]*(1-psc[:,None])
    
    s_m = Vmm - Vms[:,None]
    s_f = Vfm - Vfs[:,None]
    
    nbs = -np.inf*np.ones_like(s_m)
    
    i_pos = ((s_m > 0) & (s_f>0) )
    nbs[i_pos] = s_m[i_pos]**(gamma) * s_f[i_pos]**(1-gamma)
    
    ismar = np.any(nbs>0,axis=1)
    Vout_m, Vout_f = np.empty_like(Vms), np.empty_like(Vfs)
    
    Vout_m[~ismar] = Vms[~ismar]
    Vout_m[ismar]  = Vmm[ismar,nbs[ismar,:].argmax(axis=1)]
    Vout_f[~ismar] = Vfs[~ismar]
    Vout_f[ismar]  =  Vfm[ismar,nbs[ismar,:].argmax(axis=1)]
    
    
    return Vout_f, Vout_m#, Vms, Vfs, Vmm, Vfm, nbs
    

def EV_next(sf,trim_lvl=0.01):
    
    p_mat = np.empty((n_zf*n_zm*n_psi,zf[0].shape[0]))
    
    
   #
    
    
    for iz in range(zf[0].shape[0]):
        p_zm  = int_prob(zm[-1], mu=zf[0][iz],sig=sig_z)
        p_psi = int_prob(psi[-1],mu=0,sig=sig_psi_0)
        p_zf  = zf_mat[0][iz,:]
        #sm = sf
    
        p_vec = np.zeros(n_zf*n_zm*n_psi)
        
        
        
        for izf, p_zf_i in enumerate(p_zf):
            for izm, p_zm_i in enumerate(p_zm):
                for ipsi, p_psi_i in enumerate(p_psi):
                    p = p_zf_i*p_zm_i*p_psi_i
                    if p > trim_lvl:
                        p_vec[izf*n_zm*n_psi + izm*n_psi + ipsi] = p
        
        p_vec = p_vec / np.sum(p_vec)
        p_mat[:,iz] = p_vec
        
    
    VF_next = np.ones((1,n_zf*n_zm*n_psi))*(-1e-10)
    inds = np.where( np.any(p_mat>0,axis=1 ) )[0]
    
    for ipart in range(npart):
        sm = sf*np.exp(sig_a*v[ipart,0])
        VF_next[0,inds] = V_newmar(sf,sm,inds)[0]
        try:
            EV_f += (1/npart)*np.dot(VF_next,p_mat)
        except:
            EV_f = (1/npart)*np.dot(VF_next,p_mat)
    
    
    
    
    return EV_f#, p_mat, VF_next

print('Time elapsed is {}'.format(default_timer()-start))
pv = EV_next(2.0)
print('Time elapsed is {}'.format(default_timer()-start))
pv2 = np.array([EVnext_f(2.0,z) for z in zf[-1]]).flatten()
print('Time elapsed is {}'.format(default_timer()-start))

q  = V_newmar(np.array([0.0,0.0]),np.array([0.0,0.0]),(np.array([4,4]),np.array([4,4]),np.array([8,8])))
q2 = V_newmar(np.array([0.0,0.0]),np.array([0.0,0.0]),np.array([368,368]) )
qq = Vnext(np.array([0.0,0.0]),np.array([0.0,0.0]),zm[-1][np.array([4,4])],zf[-1][np.array([4,4])],psi[-1][np.array([8,8])])
print((q,qq))
    
    