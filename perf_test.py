#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 12:57:15 2019

@author: egorkozlov
"""

import numpy as np
from timeit import default_timer

start = default_timer()

x = np.arange(10000,dtype=np.float64).reshape(100,100)
y = np.arange(4000,dtype=np.float64)

M = x[:,:,None] - y[None,None,:]

#i_pos = M>0



#M[i_pos] = np.sqrt(M[i_pos])
#M[~i_pos] = -np.inf

'''


nrep = 100
for rep in range(nrep):
    x = np.arange(2000,dtype=np.float64)
    y = np.arange(4000,dtype=np.float64)
    
    M = x[:,None] - y[None,:]
    
    M_sq = M
    
    i_pos = M>0
    
    M_sq[i_pos] = np.sqrt(M[i_pos])
    M_sq[~i_pos] = -np.inf

stop = default_timer()
print('Time per iteration is {} sec'.format( (stop-start)/nrep) )


'''


start = default_timer()
nrep = 10
for rep in range(nrep):
    x = np.arange(10000,dtype=np.float64).reshape(100,100)
    y = np.arange(4000,dtype=np.float64)
    
    M = x[:,:,None] - y[None,None,:]
    
    M_sq = M
    
    #i_pos = M>0
    
    #M_sq[i_pos] = np.sqrt(M[i_pos])
    #M_sq[~i_pos] = -np.inf
    M = np.where(M>=0,np.sqrt(M),-np.inf)

stop = default_timer()
print('Time per iteration is {} sec'.format( (stop-start)/nrep) )


'''


from numba import jit
@jit(nopython=True)
def apply_sqrt(M):    
    M_out = np.full_like(M,-np.inf)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            for k in range(M.shape[2]):
                m_now = M[i,j,k]                 
                if m_now >= 0:
                    M_out[i,j,k] = np.sqrt(m_now)
                else:
                    break
    return M_out


M = apply_sqrt(M)



start = default_timer()
nrep = 10
for rep in range(nrep):
    x = np.arange(10000,dtype=np.float64).reshape(100,100)
    y = np.arange(4000,dtype=np.float64)
    
    M = x[:,:,None] - y[None,None,:]
    
    M = apply_sqrt(M)
    
    

stop = default_timer()
print('Time per iteration is {} sec'.format( (stop-start)/nrep) )
'''