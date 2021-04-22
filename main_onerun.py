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
    dat_moments(period=1,sampling_number=100,weighting=True,transform=2)
    
         
    #Initialize the file with parameters
    x0 = np.array([0.882993967791,0.873782221571,1.80145051864,0.367018181976,1.08466672563,0.173521741744,-0.0602557796554,1.28728416983])
    x0 = np.array([0.882993967791,0.873782221571,1.80145051864,0.367018181976,0.8466672563,0.173521741744,-0.0602557796554,1.28728416983])
    x0 = np.array([0.870537,0.454929,3.40098,0.598726,0.870482,0.123924,-0.0423595,1.04416])
    x0 = np.array([0.870537,0.454929,3.40098,0.598726,1.170482,0.123924,-0.0823595,1.04416])
       
    
    x0 = np.array([0.754453125, 0.49631347656249997, 2.281884765625,0.33953613281250006,1.042724609375, 0.1762547800292969, -0.043862304687500006,  1.1399267578124999])
    #x0 = np.array([ 0.84626953,  0.67436523,  2.03974609,  0.38600586,  1.25532227,  0.15834363, -0.07493164,  1.20098633])
    x0 = np.array([0.87737305, 0.51726074, 1.83100586,0.3410498,1.27399902,-0.09919137,-0.1061377,1.27565918])
    x0 = np.array([0.78242188,  0.61855469,  1.99375,     0.51035156,  1.11347656, -0.09506434, -0.04082031,  1.05640625])
    
    
    #Almost ok with -0.05-global+local
    x0 = np.array([0.753853,0.388988,2.20607,0.513084,1.01583,-0.0991461,-0.0832463,1.18378])
    #x0 = np.array([1.0,0.388988,2.20607,0.513084,1.01583,0.0,-0.0432463,1.18378])
    
    
    #-0.2-global only
    #x0 = np.array([0.785987,0.45995,2.02769,0.365523,1.18527,-0.133188,-0.0867234,1.38449])
    #x0 = np.array([0.779856,0.285281,1.47748,0.278688,1.20608,-0.167525,-0.0878848,1.30855])
    #x0 = np.array([0.699238,0.280692,2.25914,0.30872,1.19259,-0.133752,-0.0596778,1.16098])
    #x0 = np.array([0.65769,0.2245,1.91147,0.262793,0.999988,-0.1321,-0.0657886,1.40649])
    #x0 = np.array([0.885987,0.45995,2.02769,0.365523,1.18527,-0.133188,-0.0867234,1.38449])
    x0 = np.array([ 0.8562500000000001,0.35507812499999997,2.65390625, 0.40023437500000003, 1.0695312499999998, -0.03746101015625,  -0.056953125, 1.416796875])
   
    x0 = np.array([0.86945,0.394907,1.67524,0.362278,1.24679,-0.086061,-0.0821776,1.35918])
    x0 = np.array([0.85694,0.728984,1.99718,0.468374,1.08246,-0.076772,-0.0716345,1.19513])
    x0 = np.array([0.85694,0.688984,1.99718,0.468374,1.08246,-0.106772,-0.0716345,1.19513])


    x0 = np.array([0.838984375,0.4752685546875,1.8553222656250001,0.445654296875,1.0952392578124999,-0.10154860463867187,-0.0706201171875,1.1860839843750002])
    x0 = np.array([0.838984375,0.5652685546875,1.8553222656250001,0.445654296875,1.0952392578124999,-0.1154860463867187,-0.0706201171875,1.1860839843750002])
    x0 = np.array([0.83898437 , 0.47526855 , 1.85532227 , 0.4456543, 1.09523926 ,-0.135486, -0.07062012,  1.0928608398])
    x0 = np.array([0.83898437 , 0.47526855 , 1.85532227 , 0.4456543, 1.09523926 ,-0.135486, -0.07062012,  1.0928608398])
    #x0 = np.array([0.83898437 , 0.47526855 , 1.85532227 , 0.4456543, 1.09523926 ,-0.139486, -0.07062012,  1.0928608398])
    #x0 = np.array([0.83898437 , 0.47526855 , 1.85532227 , 0.4456543, 1.09523926 ,-0.14186, -0.07062012,  1.0928608398])


    #(199)-very very good
    x0 = np.array([0.797307,0.762312,2.17978,0.377448,1.19598,-0.151496,-0.0789312,1.14015])
    #x0 = np.array([0.977307,0.762312,2.17978,0.377448,1.19598,-0.151496,-0.0789312,1.14015])
    
    ##197 also very goo
    #x0 = np.array([0.782733,0.821696,1.98509,0.409283,1.15559,-0.158216,-0.0720787,1.13099])
   
   # x0 = np.array([0.795805,0.300991,2.43472,0.437418,1.03025,-0.115475,-0.0848882,1.24947])

  
    if system() == 'Windows':   
        path='D:/blasutto/store_model'
    else:
        path = None
    
    out, mdl, agents, res = mdl_resid(x0,return_format=['distance','models','agents','scaled residuals'],
                                      #load_from=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
                                      solve_transition=True,                                    
                                      #save_to=['mdl_save_bil.pkl','mdl_save_uni.pkl'],
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
     
     
    
        
