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
from residuals import mdl_resid
from data_moments import dat_moments
 
print('Hi!')
 
 
 
 
if __name__ == '__main__':
     
     
    #For graphs later
    graphs=True
    #Build  data moments and pickle them
    dat_moments()
    
         
    #Initialize the file with parameters
    #x0 = np.array([0.01,0.39,0.56,0.16,0.34,0.0001,0.5])
    
    #x0 = np.array([0.29936427,0.04353319,0.2627978,0.61821716,0.30178722,0.0310754,0.64637075])
    #x0 = np.array([0.31069515,0.11153252,0.26202321,0.68278782,0.28923198,0.02740209,0.61911462])
    #0.12442258, 0.01066495, 0.0364165,  0.70268823, 0.30453891
    x0 = np.array([0.25,0.75125,0.12375,0.85,0.1075,0.0075,0.25])
     
    #from p_client import compute_for_values
    #out = compute_for_values([('compute',x0)])
    #print(out)
 
    
    out, mdl = mdl_resid(x0,return_format=['distance','model'],verbose=True,calibration_report=False,draw=True,graphs=graphs)
                         
    print('Done. Residual in point x0 is {}'.format(out))
     
    
    #Indexes for the graphs
    if graphs:
        ai=30
        zfi=0
        zmi=4
        psii=3
        ti=1
        thi=10
         
        #Actual Graphs
        mdl.graph(ai,zfi,zmi,psii,ti,thi)
         
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl')
     
     
    
        
