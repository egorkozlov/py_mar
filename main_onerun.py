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
    dat_moments(period=6,sampling_number=4,transform=2)
    
         
    #Initialize the file with parameters
    
    
    
    x0 = np.array([0.701399,0.1810307,1.11501,0.643047,0.180264,0.0,0.71854,(1-0.21)/0.21])
    x0 = np.array([0.0,0.1810307,1.11501,0.543047,0.050264,0.08,0.9999,(1-0.21)/0.21])
    x0 = np.array([0.2,0.1110307,1.11501,0.543047,0.050264,0.005,-0.09])
    x0 = np.array([0.2,0.0710307,1.11501,0.543047,0.050264,0.005,-0.09])
    

             

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
                                      verbose=True,calibration_report=False,draw=graphs,graphs=graphs)
                         
    print('Done. Residual in point x0 is {}'.format(out))
     
    
    #Indexes for the graphs
    if graphs:
        ai=30
        zfi=3
        zmi=4
        psii=1
        ti=4
        thi=5
         
        #Actual Graphs
        mdl[1].graph(ai,zfi,zmi,psii,ti,thi)
        #get_ipython().magic('reset -f')
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl')
     
     
    
        
