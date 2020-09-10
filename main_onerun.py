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
     
import os
if system() != 'Darwin' and system() != 'Windows':      
    os.environ['QT_QPA_PLATFORM']='offscreen'
    
if system() == 'Darwin':
    os.environ["PATH"] += os.pathsep + '/Library/TeX/texbin/'
   

import numpy as np
from residuals import mdl_resid
from data_moments import dat_moments
 
print('Hi!')
 
 
 
 
if __name__ == '__main__':
     
     
    #import warnings
    #warnings.filterwarnings("error")
    #For graphs later
    graphs=True
    #Build  data moments and pickle them
    dat_moments(period=1,sampling_number=4,weighting=True,transform=2)
    
         
    #Initialize the file with parameters
    x0 = np.array([0.882993967791,0.873782221571,1.80145051864,0.367018181976,1.08466672563,0.173521741744,-0.0602557796554,1.28728416983])
    x0 = np.array([0.882993967791,0.873782221571,1.80145051864,0.367018181976,0.8466672563,0.173521741744,-0.0602557796554,1.28728416983])
    x0 = np.array([0.870537,0.454929,3.40098,0.598726,0.870482,0.123924,-0.0423595,1.04416])
    x0 = np.array([0.870537,0.454929,3.40098,0.598726,1.170482,0.123924,-0.0823595,1.04416])
       
    
    #Name and location of files
    if system() == 'Windows':   
        path='D:/blasutto/store_model'
    else:
        path = None
    
    out, mdl, agents, res = mdl_resid(x0,return_format=['distance','models','agents','scaled residuals'],
                                      #load_from=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      solve_transition=True,                                    
                                      save_to=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      store_path=path,
                                      verbose=True,calibration_report=False,draw=graphs,graphs=graphs,
                                      welf=False) #Switch to true for decomposition of welfare analysis
                         
    print('Done. Residual in point x0 is {}'.format(out))
     
    #assert False
    
    #Indexes for the graphs
    if graphs:
        ai=0
        zfi=2
        zmi=0
        psii=10
        ti=0
        thi=5
         
        #Actual Graphs
        mdl[0].graph(ai,zfi,zmi,psii,ti,thi)
        #get_ipython().magic('reset -f')
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl')
     
     
    
        
