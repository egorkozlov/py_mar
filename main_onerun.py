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
    #dat_moments(period=1,sampling_number=4,transform=2)
    
         
    #Initialize the file with parameters
    x0 = np.array([ 0.907796616288,0.66724274788,1.72325126224,0.257165346817,0.674319051376,0.0212082568096,-0.0641441084303,1.08579410058])
    x0 = np.array([ 0.907796616288,0.0106724274788,1.12325126224,0.257165346817,0.00000000000074319051376,0.0,0.0,1.08579410058])
    x0 = np.array([ 0.9,0.0924274788,4.12325126224,0.257165346817,0.7,-0.06,1.08579410058])
    
    #New marriage market
    x0 = np.array([ 0.833277171,0.064620668,1.90036955,0.291509255,0.510407683,0.10569942,-0.0699018671,1.19827229])
    x0 = np.array([0.696573,0.0678817,1.97245,0.353133,0.765491,0.0952097,-0.112981,1.42511])
    x0 = np.array([0.746975,0.0905286,1.89398,0.502341,0.681238,0.101228,-0.08987,1.4423])
    
    
    x0 = np.array([0.746975,0.0905286,1.89398,0.502341,1.381238,0.101228,-0.08987,1.08])


    # #Old marriage market
    x0 = np.array([0.659832,0.795317,1.76989,0.468276,0.675271,0.0964084,-0.0792931,1.13765])
    x0 = np.array([0.731934,0.705254,1.95527,0.534551,1.08686,0.0693963,-0.0864121,1.36836])
    #x0 = np.array([0.796533,0.659385,2.72502,0.467041,1.04369,0.068463,-0.127349,1.11613])
    x0 = np.array([ 0.807796616288,0.6724274788,1.72325126224,0.257165346817,1.674319051376,0.0212082568096,-0.0841441084303,1.08579410058])





    #Name and location of files
    if system() == 'Windows':   
        path='D:/blasutto/store_model'
    else:
        path = None
    
    out, mdl, agents, res = mdl_resid(x0,return_format=['distance','models','agents','scaled residuals'],
                                      #load_from=['mdl_save_bil.pkl'],#'mdl_save_uni.pkl'],
                                      solve_transition=False,                                    
                                      #save_to=['mdl_save_bil.pkl'],#'mdl_save_uni.pkl'],
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
     
     
    
        
