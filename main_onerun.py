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
     
     
    #import warnings
    #warnings.filterwarnings("error")
    #For graphs later
    graphs=True
    #Build  data moments and pickle them
    dat_moments(period=1)
    
         
    #Initialize the file with parameters
    
    
    
    x0 = np.array([0.701399,0.310307,1.11501,0.643047,0.280264,0.117317,0.21854,1.39109])


    #Name and location of files
    if system() == 'Windows':   
        path='D:/blasutto/store_model'
    else:
        path = None
    
    out, mdl, agents, res = mdl_resid(x0,return_format=['distance','models','agents','scaled residuals'],
                                      #load_from=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      #solve_transition=True,                                      
                                      #save_to=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      store_path=path,
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
        mdl[0].graph(ai,zfi,zmi,psii,ti,thi)
         
        #If you plan to use graphs only once, deselect below to save space on disk
        #os.remove('name_model.pkl')
     
     
    
        
