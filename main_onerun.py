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
import pickle
 
print('Hi!')
 
 
 
 
if __name__ == '__main__':
     
     
    #For graphs later
    graphs=True
    #Build  data moments and pickle them
    packed_stuff=dat_moments(100,weighting=False)
    with open('moments.pkl', 'wb+') as file:
        pickle.dump(packed_stuff,file)
         
    #Initialize the file with parameters
 
 
    x0 = np.array([0.04542052, 0.00727723, 0.04185624, 0.71617425, 0.28451476,0.0])#p.exp(np.array([ -1.8603,-8.1430,-1.57934,0.25130,-0.4991]))#0.08512367 -0.03874894 -0.05721577  0.57536013  0.20720013
    #0.08512367,-0.03874894,-0.05721577, 0.57536013, 0.20720013
 
 
    out, mdl = mdl_resid(x0,return_format=['distance','model'],verbose=False,calibration_report=False,draw=True,graphs=graphs)
                         
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
     
     
    
        
