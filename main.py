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
from data_moments import dat_moments
from tiktak import tiktak
import pickle
print('Hi!')

from residuals import mdl_resid

if __name__ == '__main__':
    
    
    #Build  data moments and pickle them
    packed_stuff=dat_moments(100,weighting=False)
    with open('moments.pkl', 'wb+') as file:
        pickle.dump(packed_stuff,file)
    
          
    #Create grids of parameters
    sigma_psi_g=np.linspace(0.01,0.3,1)
    sigma_psi_init_g=np.linspace(0.02,0.5,1)
    di_co_g=np.linspace(0.05,0.3,1)
    bila=np.array([False,True])
    
    
    
    #Initialize the file with parameters


    x0 = np.array([0.01,0.01,0.02,0.7,0.25,0.0001,0.5])
    lb= np.array([0.0,0.005,0.015,0.4,0.01,0.0,0.0])
    ub= np.array([1.0,1.0,0.45,1.0,0.4,0.01,1.0])
    ub[3]=min(ub[3],1.0)
    
    
    
    ##### FIRST LET'S TRY TO RUN THE FUNCTION IN FEW POINTS
    
    print('Testing the workers...')
    from p_client import compute_for_values
    pts = [lb + rs(lb.shape)*(ub-lb) for _ in range(3)]
    pts = [('compute',x) for x in pts]    
    outs = compute_for_values(pts,timeout=3600.0)
    print('Everything worked, output is {}'.format(outs))
    
    
    
    print('')
    print('')
    print('running tic tac...')
    print('')
    print('')
    
   

    #Tik Tak Optimization
    param=tiktak(200,50,10,lb,ub,mdl_resid,tole=1e-3,nelder=False,refine=False)
    
    print('f is {} and x is {}'.format(param[0],param[1]))
    
    #Now Re do the computation with graphs!
    out, mdl = mdl_resid(param[1],return_format=['distance','model'],calibration_report=False,
                         verbose=True,draw=True)
    
    
   
        

