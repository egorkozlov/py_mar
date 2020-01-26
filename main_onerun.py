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
    #dat_moments(period=2)
    
         
    #Initialize the file with parameters
    #x0 = np.array([0.01,0.39,0.56,0.16,0.34,0.0001,0.5])
    
    #x0 = np.array([0.419146,0.12785496,0.34224688,0.56194163,0.18066626,0.04687082,0.48315671])
    x0 = np.array([1.42219318, 3.52403022, 0.76691598, 0.89285592, 0.18171488,0.08080125, 0.48178463])
    
    out, mdl, agents, res = mdl_resid(x0,return_format=['distance','model','agents','scaled residuals'],
                                      #load_from=['mdl_save_uni.pkl','mdl_save_bil.pkl'],
                                      save_to='mdl_save_uni.pkl',
                                      verbose=True,calibration_report=False,draw=True,graphs=graphs)
                         
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
     
     
    
        
