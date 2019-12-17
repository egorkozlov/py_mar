#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
aCreated on Tue Sep 17 19:14:08 2019

@author: Egor Kozlov
"""





if __name__ == '__main__':
    
    try:
        from IPython import get_ipython
        get_ipython().magic('reset -f')
    except:
        pass


from platform import system
    
if system() != 'Darwin' and system() != 'Windows':   
    import os
    os.environ['QT_QPA_PLATFORM']='offscreen'


import numpy as np
from numpy.random import random_sample as rs
from tiktak import tiktak
print('Hi!')

from residuals import mdl_resid

if __name__ == '__main__':
    
    
    #If on server set Display
    
 
   
            
    #Create grids of parameters
    sigma_psi_g=np.linspace(0.01,0.3,3)
    sigma_psi_init_g=np.linspace(0.05,0.5,3)
    di_co_g=np.linspace(0.05,0.3,3)
    bila=np.array([False,True])
    
    
    
    #Initialize the file with parameters


    x0 = np.array([0.05,0.01,0.02,0.7,0.25])
    lb= np.array([0.0,0.005,0.015,0.4,0.01])
    ub= np.array([1.0,1.0,0.45,1.0,0.4])
    ub[3]=min(ub[3],1.0)
    
    
    ##### FIRST LET'S TRY TO RUN THE FUNCTION IN FEW POINTS
    
    print('Testing the workers...')
    from p_client import compute_for_values
    pts = [lb + rs(lb.shape)*(ub-lb) for _ in range(2)]
    pts = [('compute',x) for x in pts]    
    outs = compute_for_values(pts,timeout=3600.0)
    print('Everything worked, output is {}'.format(outs))
    
    
    
    x0 = np.array([0.05,0.01,0.02,0.7,0.25])#[ 0.15652793 -0.04038299  0.16475847  0.65413566  0.71109062]
    lb= x0*0.5#np.array([0.01,0.005,0.015,0.4,0.01])
    ub= x0*2.0#np.array([0.1,0.3,0.45,1.0,0.4])
    ub[3]=min(ub[3],1.0)
   
    print('')
    print('')
    print('running tic tac...')
    print('')
    print('')
    
   

    #Tik Tak Optimization
    param=tiktak(200,15,2,lb,ub,mdl_resid,tole=1e-3,nelder=False,refine=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
    
   
        

